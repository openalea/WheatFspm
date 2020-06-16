import multiprocessing as mp
import os
import time

import pandas as pd
import rearrange_graphs
import rearrange_postprocessing
import run_fspmwheat

scenarii_df = pd.read_csv(os.path.join('inputs', 'scenarii_list.csv'), index_col='Scenario')
scenarii_df['Scenario'] = scenarii_df.index
scenarii = scenarii_df.Scenario

if __name__ == '__main__':
    tstart = time.time()
    num_processes = mp.cpu_count()
    p = mp.Pool(num_processes)

    mp_solutions = p.map(run_fspmwheat.run_fspmwheat, list(scenarii))
    p.terminate()
    p.join()

    if 'Generate_Graphs' in scenarii_df.columns and any(scenarii_df.Generate_Graphs):
        rearrange_graphs.rearrange_graphs(scenarii=list(scenarii))
    if 'Run_Postprocessing' in scenarii_df.columns and any(scenarii_df.Run_Postprocessing):
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=1999, scenarii=list(scenarii))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['axes_postprocessing'], t=3499, scenarii=list(scenarii))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['organs_postprocessing'], t=3499, scenarii=list(scenarii))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['performance_indices'], scenarii=list(scenarii))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['leaf_traits'], scenarii=list(scenarii))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['canopy_kinetics_daily'], scenarii=list(scenarii))
        rearrange_postprocessing.rearrange_postprocessing(postprocessing_tables=['Conc_phloem'], scenarii=list(scenarii))

    tend = time.time()
    tmp = (tend - tstart) / 60.
    print("multiprocessing: %8.3f minutes" % tmp)
