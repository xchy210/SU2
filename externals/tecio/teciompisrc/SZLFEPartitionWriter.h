 #pragma once
#include "SZLFEZoneWriter.h"
namespace tecplot { namespace ___3933 { class ItemSetIterator; class SZLFEPartitionWriter : public SZLFEZoneWriter { public: SZLFEPartitionWriter( ItemSetIterator&                    varIter, ___4636                         zone, ___4636                         ASSERT_ONLY(___341), ___2090::___2980            ___2977, std::vector<___372> const&       ___4564, ___372                           ___4499, ___37&                         ___36, boost::shared_ptr<___1350 const> zoneInfo); virtual ~SZLFEPartitionWriter(); virtual ___2479 varMinMax(___4352 ___4336); }; }}