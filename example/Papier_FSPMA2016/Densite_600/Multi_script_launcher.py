import subprocess
import multiprocessing as mp
import time

def launch_cmd(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

def multi_script_launcher():
    commands = ['cd Asso_diffus && python main.py', 'cd Asso_direct && python main.py', 'cd Asso_mixte && python main.py',
                'cd Erect_diffus && python main.py', 'cd Erect_direct && python main.py', 'cd Erect_mixte && python main.py',
                'cd Plano_diffus && python main.py', 'cd Plano_direct && python main.py', 'cd Plano_mixte && python main.py']

    tstart = time.time()
    num_processes = mp.cpu_count()-1
    p = mp.Pool(num_processes)

    mp_solutions = p.map(launch_cmd, commands)
    p.close()
    p.join()

    tend = time.time()
    tmp = (tend - tstart) / 60.
    print("multiprocessing: %8.3f minutes" % tmp)

if __name__ == '__main__':
    multi_script_launcher()
