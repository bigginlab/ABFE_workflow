"""
This file can be used for some general standard settings.

"""

from abfe import template


std_conf = {
    #Gromacs Kernels
    "gmx_kernel_cpu": template.gmx_submit_kernels_path + "/def_gpu_job_cont.sh",
    "gmx_kernel_cpu_cont": template.gmx_submit_kernels_path + "/def_gpu_job_cont.sh",
    "gmx_kernel_gpu": template.gmx_submit_kernels_path + "/def_gpu_job_cont.sh",
    "gmx_kernel_gpu_cont": template.gmx_submit_kernels_path + "/def_gpu_job_cont.sh",


            }
