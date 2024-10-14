import pandas as pd
import argparse
import os

def merge_samples(samplesheet_path, detect_strand_path, *replicas ):
    samplesheet = pd.read_csv(samplesheet_path)
    detect_strand = pd.read_csv(detect_strand_path, sep='\t')
    groups = samplesheet['sample'].unique()
    for group in groups:
        # Get the samples for this group

        if 'rna' in samplesheet.columns:
            samples = samplesheet[samplesheet['sample'] == group]['rna'].values
        else:
            samples = samplesheet[samplesheet['sample'] == group]['fastq_1'].values

        merged_data = pd.DataFrame()

        for sample in samples:
            # Read the .tab file for this sample
            # tab_file_path = os.path.join(tab_dir, f'{sample}.tab')
            sample = os.path.basename(sample)
            sample = os.path.splitext(sample)[0]
            tab_file_path = None
            for replica in replicas:
                if sample in replica:
                    tab_file_path = replica
                    break
            
            # tab_file_path = f'{sample}.1-CIGAR.tab'
            data = pd.read_csv(tab_file_path, sep='\t', low_memory=False)
            data['srr_id'] = f'{sample}'
            
            # Check for RNA and DNA strand orientation and swap it if it is ANTI
            strand = detect_strand[ detect_strand['id'] == sample ]['strand'].values[0]
            if strand == 'ANTI':
                data['rna_strand'] = data['rna_strand'].apply(lambda x: x.replace('+', 'temp').replace('-', '+').replace('temp', '-'))
            elif strand == 'SAME':
                pass
            else:
                raise Exception("Strand orientation has not been deduced! ")
            
            if merged_data.empty:
                merged_data = data
            else:
                merged_data = pd.merge(merged_data, data, how='outer')
        merged_data.to_csv(f'{group}_merged.tab', sep='\t', index=False)
    print(replicas)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge replicas .tab files')
    parser.add_argument('samplesheet_path', type=str)
    parser.add_argument('detect_strand_result', type=str)
    parser.add_argument('tab_files', nargs='+', type=str)
    args = parser.parse_args()

    merge_samples(args.samplesheet_path, args.detect_strand_result, args.tab_files)
