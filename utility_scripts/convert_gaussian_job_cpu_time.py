# -*- coding: utf-8 -*- 
#                                                     #
#  __authors__ = Adarsh Kalikadien & Vivek Sinha      #
#  __institution__ = TU Delft                         #
#                                                     #
"""This script was used in combination with get_log_gaussian_descriptors.sh to convert the job cpu time to seconds
such that jobs' CPU time became easily comparable
"""
import sys
# convert gaussian CPU job time to seconds
days = float(sys.argv[1])
hours = float(sys.argv[2])
minutes = float(sys.argv[3])
seconds = float(sys.argv[4])

total_time = days * 24 * 3600 + hours * 3600 + minutes * 60 + seconds
print(total_time)
