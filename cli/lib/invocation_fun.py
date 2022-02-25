#!/usr/bin/env python3

def write_invocation(argv, output_path):
    import time
    import sys
    args = vars(argv).items()
    args = list(args)
    sys_cmd = sys.argv[0]

    sys_time = time.asctime(time.localtime(time.time()))

    header = "{}: {}------------------".format(sys_cmd, sys_time)

    fname = output_path + "/invocation"
    file1 = open(fname, "a")  # append mode
    file1.write(header)
    file1.write("\n")
    for arg in args:
        formatted_args = "{}: {}".format(arg[0], arg[1])
        file1.write(formatted_args)
        file1.write("\n")
    file1.write("\n")
    file1.close()