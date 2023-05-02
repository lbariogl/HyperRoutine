# o2-analysis-collision-converter -b --configuration json://configuration.json --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
# o2-analysis-track-propagation -b --configuration json://configuration.json --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
# o2-analysis-trackselection -b --configuration json://configuration.json --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
o2-analysis-event-selection -b --configuration json://configuration.json --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
o2-analysis-zdc-converter -b --configuration json://configuration.json --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
o2-analysis-lf-hypertriton-reco-task -b --configuration json://configuration.json --aod-writer-keep dangling --aod-writer-ntfmerge 300 --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000 |
o2-analysis-timestamp -b --configuration json://configuration.json --shm-segment-size 750000000000 --aod-memory-rate-limit 50000000000
