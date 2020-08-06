
  
  qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path se-33-manifest \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path demux-single-end.qza
  
  
  qiime deblur denoise-16S
  
  
  
  
## for deblur results?

  
  
 #  --type 'SampleData[Sequences]' \
 #--input-format DNAFASTAFormat 
 
qiime tools import \
  --input-path qiime_fasta  \
  --type 'FeatureData[Sequence]'\
  --output-path demux-fasta-mult.qza
  

qiime tools import \
  --input-path /u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/BLOOM.fasta  \
  --type 'FeatureData[Sequence]'\
  --output-path bloom.qza
  
  
# bloom filter

qiime tools import \
  --input-path demux-single-end.qza  \
  --type 'FeatureData[Sequence]'\
  --output-path demux_feat.qza


qiime quality-control exclude-seqs \
	--i-query-sequences demux-single-end.qza \
	--i-reference-sequences bloom.qza \
	--o-sequence-hits hits.qza \
	--o-sequence-misses misses.qza
	
## didn't work
	
# deblur

qiime quality-filter q-score \
 --i-demux demux-single-end.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza
 
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 120 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza
  
  
mv rep-seqs-deblur.qza rep-seqs.qza
mv table-deblur.qza table.qza

  
 
  
  
  
  DNAFASTAFormat