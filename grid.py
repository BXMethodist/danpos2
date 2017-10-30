import pandas as pd, numpy as np, random, os
from collections import defaultdict
from copy import deepcopy
from multiprocessing import Process, Queue
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker

def df_to_index_danpos(df, bin=3000):
    results = {}
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.split()
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def df_to_index_sk(df, bin=3000):
    results = {}
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.strip().split(',')
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             float(line[8]),
             float(line[9]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def get_stats(gene_df, df_path, criteria, bin=3000, df_function=df_to_index_danpos):
    """
    This function will select the target gene's peaks stats from the dataframe
    :param gene_df: gene:range dictionary, get the best result for each gene from transcripts range
    :param df: dataframe contain peaks from certain cutoff
    :return:
    """
    if criteria != 'skewness' and criteria != 'kurtosis':
        table_dict = df_function(df_path)
    else:
        df_function = df_to_index_sk
        table_dict = df_function(df_path)

    results = defaultdict(float)

    for k in range(gene_df.shape[0]):

        # print cur_ranges
        # if gene == 'GFI1':
        #     print gene, cur_ranges
        gene_name = gene_df.iloc[k, 0]

        chr_name, start, end = gene_df.iloc[k, 1], gene_df.iloc[k, 2], gene_df.iloc[k, 3]
        ## Here is the problem, danpos selector will consider the entire overlapped peaks
        ## The other approach is using self designed peak calling, to make sure each parameter will return different value
        cur_table = set()

        if end < start:
            mid = (start + end) / 2
            start = mid
            end = mid

        for i in range(int(start/bin), int(end/bin) + 1):
            if chr_name in table_dict and i in table_dict[chr_name]:
                table = table_dict[chr_name][i]
                cur_table = cur_table.union(table)

        if len(cur_table) == 0:
            continue

        # if gene_name == 'A2BP1':
        #     print cur_table

        selected_table = []
        for t in cur_table:
            if start < t[1] < end:
                selected_table.append(t)
            elif start < t[2] < end:
                selected_table.append(t)
            elif t[1] <= start and end <= t[2]:
                selected_table.append(t)

        # if gene_name == 'CCL13':
        #     print selected_table, start, end

        if len(selected_table) == 0:
            continue

        cur_df = pd.DataFrame(list(selected_table))

        if cur_df.shape[1] == 6:
            cur_df.columns = ['chr',
                          'start',
                          'end',
                          'width_above_cutoff',
                          'total_signal',
                          'height',]
        else:
            cur_df.columns = ['chr',
                              'start',
                              'end',
                              'width_above_cutoff',
                              'total_signal',
                              'height',
                              'skewness',
                              'kurtosis']

        if criteria == 'total_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.sum()
        elif criteria == 'height':
            cur_value = cur_df['height'].max()
        elif criteria == 'single_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.max()
        elif criteria == 'total_signal':
            cur_value = cur_df['total_signal'].sum()
        elif criteria == 'single_signal':
            cur_value = cur_df['total_signal'].max()
        #
        # # This is for kurtosis and skewness
        elif cur_df.shape[0] > 0 and criteria == 'skewness' and 'skewness' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(),'skewness']
        elif cur_df.shape[0] > 0 and criteria == 'kurtosis' and 'kurtosis' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(), 'kurtosis']

        if cur_value > results[gene_name] and criteria != 'skewness' and criteria != 'kurtosis':
            results[gene_name] = cur_value
        # this is for kurtosis and skewness
        elif criteria == 'kurtosis':
            if abs(cur_value) > abs(results[gene_name]):
                results[gene_name] = cur_value
        elif criteria == 'skewness':
            if abs(cur_value) > results[gene_name]:
                results[gene_name] = abs(cur_value)

    final = []
    for gene_name in gene_df['gene'].unique():
        final.append((gene_name, results[gene_name]))
    # print final, len(final)
    return final

def call_stats(gene_df, wigs, criteria, cutoff):
    results = defaultdict(float)

    for k in range(gene_df.shape[0]):
        gene_name = gene_df.iloc[k, 0]
        chr_name, start, end = gene_df.iloc[k, 1], gene_df.iloc[k, 2], gene_df.iloc[k, 3]
        if chr_name in wigs.genome:
            cur_signals = wigs.genome[chr_name].get_signals(start, end)
            # criterias = ['total_width']
            # criterias = ['single_width']
            # criterias = ['height']
            # criterias = ['total_signal']
            # criterias = ['single_signal']
            # criterias = ['skewness']
            # criterias = ['kurtosis']
            if criteria == 'total_width':
                cur_value = len(np.where(cur_signals > cutoff)[0])
                if results[gene_name] < cur_value:
                    results[gene_name] = cur_value
            elif criteria == 'total_width':
                pass
            elif criteria == 'single_width':
                pass
            elif criteria == 'height':
                pass
            elif criteria == 'total_signal':
                pass
            elif criteria == 'single_signal':
                pass

        else:
            continue

    final = []
    for gene_name in gene_df['gene'].unique():
        final.append((gene_name, results[gene_name]))
    return final

def random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed, number_genes=None):
    """
    random select genes for each cell type

    :param CIG_gene_table: table path
    :param exclude_list: exclude_list path
    :param all_genes: all gene GTF
    :param random_seed: random seed number to make sure each return the same random gene set
    :return:
    """
    results = []
    candidates = set()
    for gene in all_genes:
        if gene not in exclude_list and gene not in CIG_gene_df['gene']:
            candidates.add(gene)
    candidates = list(candidates)
    random.seed(random_seed)

    if number_genes is None:
        number_genes = len(candidates)

    negative_control_genes = random.sample(candidates, number_genes)

    for i in range(len(negative_control_genes)):
        try:
            results.append([negative_control_genes[i]]+list(CIG_gene_df.iloc[i, 1:]))
        except:
            results.append([negative_control_genes[i]] + list(CIG_gene_df.iloc[0, 1:]))
    negative_control_genes_df = pd.DataFrame(results)
    negative_control_genes_df.columns = CIG_gene_df.columns
    return negative_control_genes_df

def CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, criteria,
                 TSS_pos, TTS_pos, wigs):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    CIG_gene_list = list(CIG_gene_df['gene'].values)
    non_CIG_gene_list = list(non_CIG_gene_df['gene'].values)
    CIG_gene_ranges = get_range_absolute(CIG_gene_list, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos,
                                         TTS_pos)
    non_CIG_gene_ranges = get_range_absolute(non_CIG_gene_list, all_gene_GTF, up_stream_distance, down_stream_distance,
                                             TSS_pos,
                                             TTS_pos)

    CIG_results = []
    non_CIG_results = []

    for cell_type in CIG_gene_df['cell_type'].unique():
        cur_CIG_gene_list = list(CIG_gene_df[CIG_gene_df['cell_type'] == cell_type]['gene'].values)
        cur_CIG_gene_range = CIG_gene_ranges[CIG_gene_ranges['gene'].isin(cur_CIG_gene_list)]

        cur_non_CIG_gene_list = list(non_CIG_gene_df[non_CIG_gene_df['cell_type'] == cell_type]['gene'].values)
        cur_non_CIG_gene_range = non_CIG_gene_ranges[non_CIG_gene_ranges['gene'].isin(cur_non_CIG_gene_list)]

        # print cell_type, cutoff
        # print all_dfs[cell_type].keys()
        cur_df = all_dfs[cell_type][cutoff]
        # print cur_CIG_gene_range
        cur_CIG_result = get_stats(cur_CIG_gene_range, cur_df, criteria)
        cur_non_CIG_result = get_stats(cur_non_CIG_gene_range, cur_df, criteria)

        CIG_results += cur_CIG_result
        non_CIG_results += cur_non_CIG_result
    CIG_results_df = pd.DataFrame(CIG_results)
    # print CIG_results_df
    CIG_results_df.columns = ['gene', criteria]
    non_CIG_results_df = pd.DataFrame(non_CIG_results)
    non_CIG_results_df.columns = ['gene', criteria]
    return CIG_results_df, non_CIG_results_df

def CIG_selecter_all(CIG_gene_df, all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, criteria,
                     TSS_pos, TTS_pos, wigs):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    all_gene_ranges = get_range_absolute(None, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos, TTS_pos)

    all_gene_results = []

    for cell_type in CIG_gene_df['cell_type'].unique():
        cur_df = all_dfs[cell_type][cutoff]

        # overlap selecter
        cur_CIG_result = get_stats(all_gene_ranges, cur_df, criteria)

        # non overlap selecter
        # cur_CIG_result = call_stats(all_gene_ranges, wigs, criteria, cutoff)

        all_gene_results += cur_CIG_result

    all_gene_results_df = pd.DataFrame(all_gene_results)

    all_gene_results_df.columns = ['gene', criteria]

    return all_gene_results_df

"""
This module is to calculate the ranges for different combinations of parameters
Such TSS, mid of TSS and TTS. relative size for gene body
"""
def get_range_absolute(gene_list, all_gene_GTF, left_distance, right_distance, TSS_pos, TTS_pos):
    """

    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """
    if gene_list is not None:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'].isin(gene_list)]
    else:
        cur_df = all_gene_GTF

    positive_df = cur_df[cur_df['hg19.knownGene.strand'] == '+'].copy()
    negative_df = cur_df[cur_df['hg19.knownGene.strand'] == '-'].copy()

    if TSS_pos == 'TSS' and TTS_pos == 'TSS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df.loc[positive_df['left_range'] < 0, 'left_range'] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txStart'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txEnd'] - right_distance
        negative_df.loc[negative_df['left_range'] < 0, 'left_range'] = 0

    elif TSS_pos == 'TSS' and TTS_pos == 'TTS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df.loc[positive_df['left_range'] < 0, 'left_range'] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txEnd'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txStart'] - right_distance
        negative_df.loc[negative_df['left_range'] < 0, 'left_range'] = 0

    elif TSS_pos == 'MID':
        positive_df['MID'] = (positive_df['hg19.knownGene.txStart'] + positive_df['hg19.knownGene.txEnd'])/2
        negative_df['MID'] = (negative_df['hg19.knownGene.txStart'] + negative_df['hg19.knownGene.txEnd'])/2

        positive_df['left_range'] = positive_df['MID'] + left_distance
        positive_df[positive_df['left_range'] < 0] = 0
        positive_df['right_range'] = positive_df['MID'] + right_distance

        negative_df['right_range'] = negative_df['MID'] - left_distance
        negative_df['left_range'] = negative_df['MID'] - right_distance
        negative_df[negative_df['left_range'] < 0] = 0

    new_df = positive_df.append(negative_df)
    # print new_df.columns
    result_df = new_df[['hg19.kgXref.geneSymbol', 'hg19.knownGene.chrom', 'left_range', 'right_range']]
    result_df.columns = ['gene', 'chr', 'left_range', 'right_range']

    return result_df

def get_range_relative(gene_list, all_gene_GTF, left_relative, right_relative, TSS_pos, TTS_pos):
    """

    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """
    results = set()
    for gene in gene_list:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'] == gene]
        # print cur_df
        for transcript in range(cur_df.shape[0]):
            cur_chr = cur_df.iloc[transcript, 1]
            cur_strand = cur_df.iloc[transcript, 2]
            cur_start = cur_df.iloc[transcript, 3]
            cur_end = cur_df.iloc[transcript, 4]

            cur_size = cur_end - cur_start

            left_distance = cur_size * left_relative
            right_distance = cur_size * right_relative

            if TSS_pos == 'TSS' and TTS_pos == 'TSS':
                if cur_strand == '+':
                    cur_left = cur_start + left_distance
                    cur_right = cur_start + right_distance
                elif cur_strand == '-':
                    cur_right = cur_end - left_distance
                    cur_left = cur_end - right_distance
            if TSS_pos == 'TSS' and TTS_pos == 'TTS':
                if cur_strand == '+':
                    cur_left = cur_start + left_distance
                    cur_right = cur_end + right_distance
                elif cur_strand == '-':
                    cur_right = cur_end - left_distance
                    cur_left = cur_start - right_distance

            # print cur_chr, cur_left, cur_right, cur_start, cur_strand
            if cur_left > cur_right:
                continue
            results.add((gene, (cur_chr, cur_left, cur_right)))
    # print results
    results = list(results)
    result_df = pd.DataFrame(results)
    result_df.columns = ['gene', 'range']
    return result_df

def next_grid(all_range, range_step, cur_step, cur_center, reduction_step, step_limit, search):
    new_step = cur_step/reduction_step
    if new_step < step_limit:
        new_step = step_limit
    center_index = all_range.index(cur_center)

    # print center_index
    if search:
        if range_step == cur_step:
            start_index = center_index - 2 if center_index - 2 >= 0 else 0
            end_index = center_index + 3 if center_index + 3 <= len(all_range) else len(all_range)
            # print start_index, end_index
            return all_range[start_index:end_index], new_step
        else:
            start_index = center_index - int(cur_step*2/range_step) if center_index - int(cur_step*2/range_step) >=0 else 0
            end_index = center_index + int(cur_step*2/range_step) + 1 if center_index + int(cur_step*2/range_step) + 1 <= len(all_range) else len(all_range)
            # print start_index, end_index

            return all_range[start_index:end_index][::int(new_step/range_step)], new_step
    else:
        start_index = center_index - 4 if center_index - 4 >= 0 else 0
        end_index = center_index + 5 if center_index + 5 <= len(all_range) else len(all_range)
        return all_range[start_index:end_index][::int(new_step / range_step)], new_step

def grid_search(CIG_gene_df, non_CIG_gene_df,
                all_gene_GTF, all_dfs, criteria, marker, cost_function,
                TSS_pos, TTS_pos,
                up_stream_distance_range, window_size_range, cutoff_range,
                up_stream_distance_grid=10000, window_size_grid=10000, cutoff_grid=10,
                up_stream_distance_range_step=1000, window_size_range_step=1000, cutoff_range_step=1,
                up_stream_distance_step=2, window_size_step=2, cutoff_step=2,
                up_stream_distance_limit=1000, window_size_limit=1000, cutoff_limit=1,
                process=8, wigs=None):

    ## track the parameters and logP path
    path = []

    new_up_stream_distance_range = deepcopy(up_stream_distance_range)
    new_window_size_range = deepcopy(window_size_range)
    new_cutoff_range = deepcopy(cutoff_range)

    search = True

    # iter = 1
    # count = 0

    grid_up_stream_distance_range = new_up_stream_distance_range[
                                    up_stream_distance_grid / up_stream_distance_range_step / 2::(
                                    up_stream_distance_grid / up_stream_distance_range_step)]
    grid_window_size_range = new_window_size_range[
                             window_size_grid / window_size_range_step / 2::(window_size_grid / window_size_range_step)]
    grid_cutoff_range = new_cutoff_range[int(cutoff_grid / cutoff_range_step / 2)::int(cutoff_grid / cutoff_range_step)]

    # grid_up_stream_distance_range = [4000]
    # grid_window_size_range = range(1000, 500000)
    # grid_cutoff_range = [10]

    past_path = {}

    while True:
        # determine the new grid

        print grid_up_stream_distance_range
        print grid_window_size_range
        print grid_cutoff_range

        # if count >= iter:
        #     break

        combinations = np.stack(np.meshgrid(grid_up_stream_distance_range, grid_window_size_range, grid_cutoff_range), -1).reshape(-1, 3)

        # print combinations
        if TSS_pos == TTS_pos:
            new_combinations = []
            for comb in combinations:
                if comb[1] > comb[0]:
                    new_combinations.append(comb)
            combinations = np.asarray(new_combinations)
        print combinations.shape, 'current grid size is '

        # combinations = np.asarray([[2000, 3000, 8]])

        # print combinations
        #
        # return

        best_log_P = None
        best_comb = None

        chunks = []
        cur_index = 0
        reminder = len(combinations) % process
        chunk_size = len(combinations) / process
        for i in range(process):
            if reminder > 0:
                chunks.append(combinations[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
                cur_index += 1
                reminder -= 1
            else:
                chunks.append(combinations[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

        total_chunk_size = 0
        for chunk in chunks:
            total_chunk_size += len(chunk)
        if total_chunk_size != len(combinations):
            print 'multiple processes chunk size is not correct'
            return None

        queue = Queue()
        processes = []

        cur_past_path = deepcopy(past_path)
        for i in range(process):
            cur_chunk = chunks[i]
            p = Process(target=CIG_process,
                        args=(queue, cur_chunk, CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                              all_dfs, criteria, cur_past_path, cost_function, marker, TSS_pos, TTS_pos,
                              wigs))
            processes.append(p)
            p.start()

        cur_path = []

        for i in range(process):
            cur_best_logP, cur_best_comb, cur_process_path = queue.get()
            cur_path += cur_process_path
            # print cur_best_logP, cur_best_logP< best_log_P, best_log_P is None or cur_best_logP < best_log_P
            if (best_log_P is None or cur_best_logP < best_log_P) and cur_best_logP is not None:
                best_log_P = cur_best_logP
                best_comb = cur_best_comb
            print best_comb, best_log_P
        print best_comb, best_log_P
        for p in processes:
            p.join()

        for trace in cur_path:
            if tuple(trace[0]) not in past_path:
                past_path[tuple(trace[0])] = trace[1]

        path += cur_path

        if not search:
            break

        # break
        if up_stream_distance_grid == up_stream_distance_limit and window_size_grid == window_size_limit and cutoff_grid == cutoff_limit:
            search = False
            # break
        # use best comb to determine new range

        up_stream_center, window_size_center, cutoff_center = best_comb
        grid_up_stream_distance_range, up_stream_distance_grid = next_grid(up_stream_distance_range,
                                                                           up_stream_distance_limit,
                                                                           up_stream_distance_grid,
                                                                           up_stream_center,
                                                                           up_stream_distance_step,
                                                                           up_stream_distance_limit,
                                                                           search)
        grid_window_size_range, window_size_grid = next_grid(window_size_range,
                                           window_size_limit,
                                           window_size_grid,
                                           window_size_center,
                                           window_size_step,
                                           window_size_limit,
                                           search)
        grid_cutoff_range, cutoff_grid = next_grid(cutoff_range,
                                      cutoff_limit,
                                      cutoff_grid,
                                      cutoff_center,
                                      cutoff_step,
                                      cutoff_limit,
                                      search)

        # break

    return path


def CIG_process(queue, combinations, CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                all_dfs, criteria, cur_past_path, cost_function, marker, TSS_pos, TTS_pos,
                wigs):
    path = []
    best_log_P = None
    best_comb = None

    for comb in combinations:
        # print comb[0], comb[1], comb[1] <= comb[0]

        cur_up_stream_distance, cur_window_size, cur_cutoff = comb
        if (cur_up_stream_distance, cur_window_size, cur_cutoff) in cur_past_path:
            cur_logP = cur_past_path[(cur_up_stream_distance, cur_window_size, cur_cutoff)]
        else:
            cur_logP = cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                     cur_up_stream_distance, cur_window_size,
                                     all_dfs, cur_cutoff, criteria, marker, TSS_pos, TTS_pos,
                                     wigs)
        print comb, cur_logP
        path.append((comb, cur_logP))
        if (best_log_P is None or best_log_P > cur_logP) and cur_logP is not None:
            best_log_P = cur_logP
            best_comb = comb
    # print best_comb, best_log_P, 'current process best value is'
    #     break
    # print 'for current process', best_log_P, best_comb
    queue.put((best_log_P, best_comb, path))
    return

def logP_wilcoxon(groupA, groupB, bags=1):
    """
    return the negative log P value for two groups
    :param groupA:
    :param groupB:
    :return:
    """
    p_values = []

    for i in range(bags):
        # cur_groupA = np.random.choice(groupA, int(len(groupA)*0.75), replace=True)
        # cur_groupB = np.random.choice(groupB, int(len(groupB)*0.75), replace=True)
        cur_groupA = groupA
        cur_groupB = groupB
        try:
            rank_diff, p = stats.mannwhitneyu(cur_groupA, cur_groupB, alternative='less')
            # print p
            # There is a problem in multiple processing. for some reason, the function will return a very small number
            # However, it could pass the test without problem if just run the corresponding dataframe
            # if np.log10(p) < -80:
            #     p_values.append(0)

            p_values.append(np.log10(p))
        except:
            p_values.append(0)

    return np.mean(p_values)

def logP_fisher(gene_df, all_stat_df, criteria, top_enrich=500, ascending=False):
    total_genes = all_stat_df.shape[0]
    sort_result = all_stat_df.sort_values(by=[criteria], ascending=ascending)
    top_genes = sort_result['gene'].tolist()[:top_enrich]
    overlap = len(set(top_genes).intersection(set(gene_df['gene'].tolist())))
    not_overlap = top_enrich - overlap
    # print overlap, top_enrich, not_overlap, total_genes
    p = stats.fisher_exact([[overlap, not_overlap], [not_overlap, total_genes - 2 * top_enrich + overlap]],
                           alternative='greater')[1]
    return np.log10(p)

def wilcoxon_cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                          cur_up_stream_distance, cur_down_stream_distance,
                          all_dfs, cur_cutoff, criteria, marker,
                          TSS_pos, TTS_pos, wigs):
    cur_CIG_results_df, cur_non_CIG_results_df = CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                                              cur_up_stream_distance, cur_down_stream_distance,
                                                              all_dfs, cur_cutoff, criteria,
                                                              TSS_pos, TTS_pos,
                                                              wigs)

    # print 'wilcoxon'
    print cur_CIG_results_df[criteria].mean(),  cur_non_CIG_results_df[criteria].mean(), 'average'

    if cur_CIG_results_df[criteria].mean() < cur_non_CIG_results_df[criteria].mean():
        cur_logP = logP_wilcoxon(cur_CIG_results_df[criteria],
                                 cur_non_CIG_results_df[criteria])
    else:
        cur_logP = logP_wilcoxon(cur_non_CIG_results_df[criteria],
                                 cur_CIG_results_df[criteria])

    return cur_logP

def fisher_cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                         cur_up_stream_distance, cur_down_stream_distance,
                         all_dfs, cur_cutoff, criteria, marker,
                         TSS_pos, TTS_pos, wigs):
    all_gene_results_df = CIG_selecter_all(CIG_gene_df, all_gene_GTF, cur_up_stream_distance, cur_down_stream_distance,
                                           all_dfs, cur_cutoff, criteria,
                                           TSS_pos, TTS_pos, wigs)
    # print 'fisher'
    cur_logP = logP_fisher(CIG_gene_df, all_gene_results_df, criteria, top_enrich=500)

    return cur_logP

def reformat(directory='./csv/'):
    paths = [x for x in os.listdir(directory) if x.find('change') == -1 and x.endswith('.csv')]
    for path in paths:
        results = []
        print path
        df = pd.read_csv(directory+ path, index_col=0)

        print df

        df.columns = ['p','logP']

        for i in range(df.shape[0]):
           info = df.ix[i, 'p']
           info = info.replace('[', '')
           info = info.replace(']', '')
           cur_result = info.split() + [df.ix[i, 'logP']]
           results.append(cur_result)

        df = pd.DataFrame(results)
        df.columns = ['upstream', 'downstream', 'height', 'logP']
        # df = df[df['logP'] > -82]
        df.to_csv(directory+path, index=None)

def get_best(directory='./csv/'):
    parameters = {'upstream': range(-1000000, 1000000, 1000), 'downstream': range(-1000000, 1000000, 1000), 'height': [0,0.25, 0.5,0.75, 1.0,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 ,4.5, 5.0] + range(5, 301)}
    # parameters = {'upstream': range(-1000000, 1000000, 1000), 'downstream': range(-1000000, 1000000, 1000),
    #               'height': range(5, 301)}
    markers = ['h3k4me1', 'h3k27ac', 'h3k4me3', 'h3k27me3']
    # markers = ['h3k4me1', 'h3k4me3', 'h3k27ac']
    # markers = ['h3k27me3']
    features = ['total_width', 'height', 'total_signal', 'kurtosis', 'skewness']
#['total_width', 'single_width',  'total_signal', 'kurtosis', 'skewness', 'single_signal',]
    for marker in markers:
        for feature in features:
            for parameter in parameters.keys():
                results = []
                df = pd.read_csv(directory + 'grid_path_'+marker+'_'+feature+'.csv')
                for p in parameters[parameter]:
                    cur_df = df[df[parameter] == p]
                    for other_parameter in parameters.keys():
                        if other_parameter != parameter:
                            # print parameters[other_parameter]
                            cur_df = cur_df[cur_df[other_parameter].isin(parameters[other_parameter])]
                    if cur_df.shape[0] == 0:
                        continue
                    results.append((p, cur_df['logP'].min()*-1))
                df = pd.DataFrame(results)
                df.columns = [parameter, 'best_logP']
                df.to_csv(directory+marker+'_'+feature+'_'+parameter+'.csv', index=None)

def get_best_parameter(directory='./csv/'):
    results = {}
    paths = [x for x in os.listdir(directory) if x.find('change') == -1 and x.endswith('.csv')]
    for path in paths:
        df = pd.read_csv(directory+path)
        marker = path.split('_')[2]
        info = path[:-4].split('_')
        feature = '_'.join(info[3:])
        best_p = df['logP'].min()
        cur_df = df[df['logP'] == best_p]
        if (marker, feature) in results:
            if best_p < results[(marker, feature)][-1]:
                results[(marker, feature)] = [tuple(x) for x in cur_df.to_records(index=False)][0]
        else:
            results[(marker, feature)] = [tuple(x) for x in cur_df.to_records(index=False)][0]
    final = []
    for key, value in results.items():
        final.append([key[0], key[1]] + list(results[key]))

    final_df = pd.DataFrame(final)
    final_df.columns = ['marker', 'feature', 'upstream', 'downstream', 'height', 'logP']
    final_df.to_csv('best_parameters_CIG.csv', index=None)

def generate_plot(df_path, marker, verbose=True):
    df = pd.read_excel(df_path)
    df = df[df[marker]!=0]
    if marker == 'upstream':
        df[marker] = df[marker] * -1
    df = df.set_index([marker])
    columns = df.columns
    colors = ['red', 'blue', 'green', 'purple', 'yellow']
    shapes = ['o', 'D','^', 'x']
    labels = []
    handals = []

    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    if marker != 'height':
        ax.set_xscale('symlog')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

    for i in range(0, len(columns)):
        shape = shapes[i]
        # ax = df.plot.scatter(columns[0], columns[i], color=colors[i], xlim=(df[columns[0]].min(), df[columns[0]].max()),
        #                      ylim=(0, df[columns[i]].max()),
        #                       edgecolor='', ax=ax)
        if marker == 'downstream':
            ax = plt.scatter(df.index/1000, df[columns[i]], color=colors[i], edgecolor='', marker=shape)
            mask = np.isfinite(df[columns[i]])
            plt.plot(df.index[mask]/1000, df[columns[i]][mask], color=colors[i], linestyle='-')

        elif marker == 'upstream':
            ax = plt.scatter(df.index / 1000, df[columns[i]], color=colors[i], edgecolor='', marker=shape)
            mask = np.isfinite(df[columns[i]])
            plt.plot(df.index[mask]/1000, df[columns[i]][mask], color=colors[i], linestyle='-')
        else:
            ax = plt.scatter(df.index, df[columns[i]], color=colors[i], edgecolor='', marker=shape)
            mask = np.isfinite(df[columns[i]])
            plt.plot(df.index[mask], df[columns[i]][mask], color=colors[i], linestyle='-')

        handals.append(ax)
        labels.append(columns[i])

    # print handles, labels
    # ax.legend(handles, labels, loc='best')
    font = {'fontname': 'Helvetica','fontsize':10}

    if marker == 'downstream':
        plt.xlabel('Downstream distance (kb)', **font)
        plt.xlim((df.index.min()/1000, df.index.max()/1000))
    elif marker == 'upstream':
        plt.xlabel('Upstream distance (kb)', **font)
        plt.xlim(df.index.max() / 1000, df.index.min() / 1000)
    else:
        plt.xlabel('Height cutoff', **font)
        # plt.xlim((np.log10(df.index.min()), np.log10(df.index.max())))
        # plt.xlim((0, df.index.max()))
        plt.xlim((0, 100))
    #
    plt.ylabel('-log10 enrich P', **font)
    # print labels
    # print ((0, df.max().max()))
    plt.ylim(((0, (int(df.max().max())/10+1)*10)))
    plt.xticks(**font)
    plt.yticks(**font)
    print labels
    if verbose:
        plt.legend(handals, labels, scatterpoints=1, loc=9, fontsize=10, ncol=2)

    plt.savefig(df_path.replace('.xlsx', '.pdf'))
    plt.close('all')

def group_plot(marker, ranges, path):
    csvs = [x for x in os.listdir(path) if x.find(marker)!= -1 and not x.startswith('grid') and not x.endswith('.xlsx') and x.endswith('.csv')]
    # print csvs
    csvs = [x for x in csvs if x.split('_')[2].find(marker)!=-1]
    # print csvs, marker
    dfs = [pd.read_csv(path+'/'+ c) for c in csvs]
    # names = [c.split('_')[0]+'_'+''.join(c.split('_')[2:]).replace('onco', 'oncogene').replace('sup', 'suppressor').replace('inputgenebody', '_genebody').replace('input', "_TSS").replace(marker, '').replace('.csv','').replace('width','') for c in csvs]
    names = [c.split('_')[0] for c in csvs]
    print names
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop_duplicates()
        dfs[i] = dfs[i].set_index([marker])
    results = []
    for r in ranges:
        cur_results = [r]
        for i in range(len(dfs)):
            try:
                cur_results.append(dfs[i].ix[r, 'best_logP'])
            except:
                cur_results.append(None)
        results.append(cur_results)
    result_df = pd.DataFrame(results)
    result_df.columns=[marker] + names
    result_df.to_excel(path+'/'+marker+'.xlsx', index=False)

    generate_plot(path+'/'+marker+'.xlsx', marker)
    return

def grid(target_table1, target_table2,
         GTF, all_dfs, feature, cost_function,
         up_stream_range, up_stream_grid, up_stream_distance_range_step, up_stream_grid_step, up_stream_grid_limit,
         down_stream_range, down_stream_grid, down_stream_distance_range_step, down_stream_grid_step, down_stream_grid_limit,
         height_range, height_grid, height_range_step, height_grid_step, height_grid_limit,
         TSS_pos, TTS_pos, output_prefix, output_path,
         process=8):
    if cost_function == 'fisher':
        grid_search(target_table1, target_table2,
                    GTF, all_dfs,  criteria=feature, marker='', cost_function=fisher_cost_function,
                    TSS_pos=TSS_pos, TTS_pos=TTS_pos,
                    up_stream_distance_range=up_stream_range, window_size_range=down_stream_range, cutoff_range=height_range,
                    up_stream_distance_grid=up_stream_grid, window_size_grid=down_stream_grid, cutoff_grid=height_grid,
                    up_stream_distance_range_step=up_stream_distance_range_step, window_size_range_step=down_stream_distance_range_step, cutoff_range_step=height_range_step,
                    up_stream_distance_step=up_stream_grid_step, window_size_step=down_stream_grid_step, cutoff_step=height_grid_step,
                    up_stream_distance_limit=up_stream_grid_limit, window_size_limit=down_stream_grid_limit, cutoff_limit=height_grid_limit,
                    process=process)
    elif cost_function == 'wilcoxon':
        grid_search(target_table1, target_table2,
                    GTF, all_dfs, criteria=feature, marker='', cost_function=wilcoxon_cost_function,
                    TSS_pos=TSS_pos, TTS_pos=TTS_pos,
                    up_stream_distance_range=up_stream_range, window_size_range=down_stream_range,
                    cutoff_range=height_range,
                    up_stream_distance_grid=up_stream_grid, window_size_grid=down_stream_grid, cutoff_grid=height_grid,
                    up_stream_distance_range_step=up_stream_distance_range_step,
                    window_size_range_step=down_stream_distance_range_step, cutoff_range_step=height_range_step,
                    up_stream_distance_step=up_stream_grid_step, window_size_step=down_stream_grid_step,
                    cutoff_step=height_grid_step,
                    up_stream_distance_limit=up_stream_grid_limit, window_size_limit=down_stream_grid_limit,
                    cutoff_limit=height_grid_limit,
                    process=process)
