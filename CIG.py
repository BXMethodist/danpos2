### This is the pipeline for CIG project

import os, pandas as pd

# Step 1: Run Danpos

def danpos_input_ENC(sample_id, input_id, node_id):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + sample_id + '.bowtie' +' -b '+ input_id +".bowtie" + danpos_parameters

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    pbs.write(cmd + '\n')
    pbs.close()
    # os.system('qsub ' + sample_id + ".pbs")
    return

def danpos_inputs(sample_id):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + sample_id + '.bowtie' +' -b '+ os.getcwd() + '/' + sample_id+'_input' + danpos_parameters

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    pbs.write(cmd + '\n')
    pbs.close()
    # os.system('qsub ' + sample_id + ".pbs")
    return

def RunDanposENC(sample_input="ENC_H3K4me3_sample_pairs.csv"):
    """
    :param sample_input: file path for sample input pair
    :return:
    """
    sample_input_df = pd.read_csv(sample_input, index_col=None)

    finished = [x for x in os.listdir('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/') if
                x.endswith('bowtie')]

    ### This is the special part to re-run error bowtie file
    incomplete_f = open('incomplete.txt', 'r')

    incomplete =[line.strip() for line in incomplete_f]
    incomplete_f.close()
    incomplete_set = set()
    ###

    names = {}
    for i in range(len(finished)):
        x = finished[i]
        x = x[:-7]
        if x.find('_') != -1:
            new_name = '_'.join(tuple(sorted(x.split('_'))))
        else:
            new_name = x
        names[new_name] = x

    print names

    nodes = [1, 2, 3, 4, 5, 6]

    count = 0
    for i in range(sample_input_df.shape[0])[0:50]:
        node_id = nodes[i%6]
        sample_id = sample_input_df.ix[i, 'sample']
        input_id = sample_input_df.ix[i, 'input']

    ### This is the special part to re-run error bowtie file
        # if sample_id in incomplete or input_id in incomplete:
        #     incomplete_set.add(sample_id)
        #     incomplete_set.add(input_id)
        #     if sample_id in incomplete:
        #         print sample_id, 'incomplete'
        #     if input_id in incomplete:
        #         print input_id, 'incomplete'
        #     # os.system('rm -r '+sample_id.strip()+'/')
        #     pass
        #     # continue
        # else:
        #     continue
    # print len(incomplete_set)
    ###

    # '''
        if pd.isnull(input_id):
            continue

        if sample_id.find('_') != -1:
            if '_'.join(tuple(sorted(sample_id.split('_')))) in names:
                sample_id = names['_'.join(tuple(sorted(sample_id.split('_'))))]
            else:
                print sample_id
                pass

        if input_id.find('_') != -1 and input_id.find(';') == -1:
            if '_'.join(tuple(sorted(input_id.split('_')))) in names:
                input_id = names['_'.join(tuple(sorted(input_id.split('_'))))]
            else:
                print input_id
                pass


        if input_id.find(';') != -1:
            inputs = [x.strip() for x in input_id.split(';')]
            if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + sample_id + '.bowtie') and \
                    (not os.path.isfile('./' + sample_id + '.bowtie')):
                os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + sample_id + '.bowtie ' + sample_id + '.bowtie')
                pass
            else:
                print sample_id, 'sample bowtie not found'
                pass
            # os.system('mkdir '+sample_id)
            for input_file in inputs:
                if input_file.find('_') != -1:
                    input_file = names['_'.join(tuple(sorted(input_file.split('_'))))]
                if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+input_file+'.bowtie') and \
                        (not os.path.isfile('./' + input_file + '.bowtie')):
                    os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + input_file + '.bowtie ' + '/' + sample_id+'_input/'+ input_file + '.bowtie')
                    pass
                else:
                    print input_file, 'input bowtie not found'

            danpos_inputs(sample_id)
            count += 1

        elif os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+sample_id+'.bowtie') \
            and os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+input_id+'.bowtie'):

            if (not os.path.isfile('./' + sample_id + '.bowtie')):
                os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/"+sample_id+'.bowtie '+sample_id+'.bowtie')
            if (not os.path.isfile('./' + input_id + '.bowtie')):
                os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + input_id + '.bowtie ' + input_id + '.bowtie')
            danpos_input_ENC(sample_id, input_id, node_id)
            count +=1
            pass
        else:
            print sample_id, input_id, 'sample or input bowtie not found'
            continue

    for pbs in [x for x in os.listdir('./') if x.endswith('.pbs')]:
        os.system('qsub '+pbs)
    print count

# Step 2: Run Selector

# python /archive/tmhkxc48/tools/danposTemp/danpos.py selector /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.regions.xls --genicSelector TSS:-3000:3000 --gene_file notch.txt_profilelist.txt --gene_out ./selector/notch.xlsx --out ./selector/notch_pathway.xlsx


# Step 3: Calculate the total width, max height..total signal...for each gene



# Step 4: Calculate the negative log P
