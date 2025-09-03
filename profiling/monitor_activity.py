import matplotlib.pyplot as plt
import numpy as np 
import sys
import logging
import os
import yaml

logger = logging.getLogger(__name__)

def main():
    test_name = str(sys.argv[1])

    # set up paths        
    config_file = '../config/config.yml'
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)
    else:
        warnings.warn(f"! Configuration file '{config_file}' not found. Using default settings.", UserWarning)
        config = {}

    default_result_path = f"../results/{test_name}/profiling/"
    result_path = config.get("result_path")
    if result_path is None:
        warnings.warn(f"! 'result_path' not found in '{config_file}'. Using default: '{default_result_path}'.", UserWarning)
        result_path = default_result_path
    else:
        result_path += f"{test_name}/profiling/"
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    activity_file = f'{result_path}/activity.txt'
    logging_file = f'{result_path}{test_name}.log'
    logging.basicConfig(
        filename=logging_file,
        level=logging.INFO,
        filemode='a+',
        format="%(message)s",
    )
    logger.info('\nReading profiling results')
    
    time, cpu_perc, real_mem, virtual_mem, read_count, write_count, read_bytes, write_bytes = np.loadtxt(activity_file, unpack=True)
    
    # Print stats
    
    logger.info(f'Process took {time.max():.2e} seconds.')
    logger.info('\nCPU:')
    logger.info(f' Max: {cpu_perc.max():.2f}% ({cpu_perc.max()/100.:.1f} threads).')
    logger.info(f' Average: {cpu_perc.mean():.2f}% ({cpu_perc.mean()/100.:.1f} threads).')
    logger.info('\nRAM:')
    # logger.info(f'Max real (virtual) memory used: {real_mem.max():.2e}MB ({virtual_mem.max():.2e}MB)')
    # logger.info(f'... per CPU: {real_mem.max()/cpu_perc.max()/100.:.2f}MB ({virtual_mem.max()/cpu_perc.max()/100.:.2f}MB)')
    logger.info(f' Max: {real_mem.max():.2e}MB')
    logger.info(f' ... per CPU: {real_mem.max()/cpu_perc.max()/100.:.2f}MB')
    logger.info(f' Average: {real_mem.mean():.2e}MB')
    logger.info('\nI/O:')
    logger.info(' Read:')
    logger.info(f'  Max count: {read_count.max():.2e} ({read_bytes.max()/1e6:.1e}MB)')
    logger.info(f'  Average count: {read_count.mean():.2e} ({read_bytes.mean()/1e6:.1e}MB)')
    logger.info(' Write:')
    logger.info(f'  Max count: {write_count.max():.2e} ({write_bytes.max()/1e6:.1e}MB)')
    logger.info(f'  Average count: {write_count.mean():.2e} ({write_bytes.mean()/1e6:.1e}MB)')
    
    # Plot stats

    # CPU
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    ax.plot(time, cpu_perc, "-", lw=1, color="C0")
    ax.set_ylabel("CPU (%)", color="C0")
    ax.set_xlabel("time (s)")
    ax.set_ylim(0.0, max(cpu_perc) * 1.2)
    ax.grid()
    fig.tight_layout()
    fig.savefig(activity_file[:-4]+'_cpu.png', dpi=220)

    # RAM
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    ax.plot(time, real_mem, "-", lw=1, color="C0", label='Real (RAM)')
    ax.plot(time, virtual_mem, "--", lw=1, color="C0", label='Virtual')
    ax.set_xlabel("time (s)")
    ax.set_ylim(0.0, max(np.r_[real_mem, virtual_mem]) * 1.2)
    ax.set_ylabel("Memory (MB)")
    ax.legend(loc='best')
    ax.grid()
    fig.tight_layout()
    fig.savefig(activity_file[:-4]+'_ram.png', dpi=220)

    # IO
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 4))
    ax.plot(time, read_count, color='C0')
    ax.plot(time, write_count, color='C0', ls='--')
    ax.set_ylabel("Count", color="C0")
    ax.set_xlabel("time (s)")
    
    ax1 = ax.twinx()
    ax1.plot(time, read_bytes/1e6, color='C1')
    ax1.plot(time, write_bytes/1e6, color='C1', ls='--')
    # ax1.set_ylim(0.0, max(np.r_[read_bytes, write_bytes]) * 1.2)
    ax1.set_ylabel("MB", color="C1")
    ax1.plot([], [], color='k', label='Read')
    ax1.plot([], [], color='k', label='Write', ls='--')
    ax1.legend(loc='best')
    
    ax.grid()
    fig.tight_layout()
    fig.savefig(activity_file[:-4]+'_io.png', dpi=220)

    logger.info('Done.')

if __name__ == '__main__':
    main()