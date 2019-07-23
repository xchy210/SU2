import shutil
from enum import IntEnum

import torch
import pysu2ad as pysu2
from mpi4py import MPI


class RunCode(IntEnum):
    # Run codes for worker processes
    STOP = -1
    RUN_FORWARD = 0
    RUN_ADJOINT = 1


def run_forward(comm, forward_driver, inputs, num_diff_outputs):
    # TODO Add checks to make sure these are run before Preprocess?
    for i, x in enumerate(inputs):
        forward_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
    forward_driver.ApplyDiff_Inputs_Vars()

    forward_driver.Preprocess(0)
    forward_driver.Run()
    forward_driver.Output(0)
    comm.Barrier()
    # TODO Way to get results in-memory, without writing to file?
    if comm.Get_rank() == 0:
        shutil.move("./restart_flow.dat", "./solution_flow.dat")

    outputs = [torch.tensor(forward_driver.GetDiff_Outputs_Vars(i))
               for i in range(num_diff_outputs)]

    for i in range(num_diff_outputs):
        if outputs[i].shape[0] > 1:
            # if dealing with full-grid, reorder according to GlobalIndex
            if comm.Get_size() > 1:
                # gather outputs in rank 0 if more than one rank
                outputs[i] = comm.gather(outputs[i], root=0)
                global_inds = comm.gather(forward_driver.GetAllGlobalIndices(), root=0)
                if comm.Get_rank() == 0:
                    outputs[i] = torch.cat(outputs[i])
                    global_inds = list(sum(global_inds, tuple()))  # join tuples
            else:
                global_inds = list(forward_driver.GetAllGlobalIndices())

            if comm.Get_rank() == 0:
                # TODO Make the list integers on the C side
                global_inds = torch.tensor(global_inds, dtype=torch.long)
                assert outputs[i].shape[0] == len(global_inds), \
                    'Only full grid outputs supported by now (besides scalars).'
                outputs[i][global_inds] = outputs[i].clone()  # order by global_inds
            else:
                outputs[i] = None

    forward_driver.Postprocessing()
    return tuple(outputs)


def run_adjoint(comm, adjoint_driver, inputs, grad_outputs):
    # TODO Add checks to make sure these are run before Preprocess
    for i, x in enumerate(inputs):
        adjoint_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
    adjoint_driver.ApplyDiff_Inputs_Vars()
    for i, g in enumerate(grad_outputs):
        adjoint_driver.SetBackprop_Derivs(g.flatten().tolist(), i)

    adjoint_driver.Preprocess(0)
    adjoint_driver.Run()
    comm.Barrier()
    grads = tuple(torch.tensor(adjoint_driver.GetTotal_Sens_Diff_Inputs(i))
                  for i in range(adjoint_driver.GetnDiff_Inputs()))
    adjoint_driver.Postprocessing()
    return grads


def main():
    intercomm = MPI.Comm.Get_parent()
    comm = intercomm.Merge(high=True)

    num_zones, dims, num_diff_outputs = comm.bcast(None, root=0)

    forward_driver, adjoint_driver = None, None
    inputs = None
    while True:
        run_type = comm.bcast(None, root=0)
        if run_type == RunCode.STOP:
            break

        config = comm.bcast(None, root=0)

        if run_type == RunCode.RUN_FORWARD:
            inputs = comm.bcast(None, root=0)
            forward_driver = pysu2.CSinglezoneDriver(config, num_zones, dims, comm)
            run_forward(comm, forward_driver, inputs, num_diff_outputs)

        elif run_type == RunCode.RUN_ADJOINT:
            assert inputs is not None, 'Run forward simulation before running the adjoint.'
            grad_outputs = comm.bcast(None, root=0)
            adjoint_driver = pysu2.CDiscAdjSinglezoneDriver(config, num_zones, dims, comm)
            run_adjoint(comm, adjoint_driver, inputs, grad_outputs)

    # TODO Disconnects hanging on cluster
    # comm.Disconnect()
    # intercomm.Disconnect()


if __name__ == '__main__':
    # will run this when run with mpirun
    main()