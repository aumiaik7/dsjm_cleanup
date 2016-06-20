#!/usr/bin/env python
import sys
import os
import datetime
import shutil
import subprocess
import gcolor
from Gcolor_OptionParser import Gcolor_OptionParser as myopt

def getCommand(exePath,fileName):
    command = exePath + " -i " +  fileName

    return command

def main():
    print "Hello World"

    parser = myopt()
    parser.gcolor_init()

    (options,args) = parser.parse_args()

    exePath = "./ColumnsAppearSorted"


    turn = 1

    for fileName in gcolor.getmtxfiles(options.root):
        print fileName
        command = getCommand(exePath,fileName)
        print command
        os.system(command)
        os.system("pkill -9 ColumnsAppearSorted")

if __name__ == "__main__":
    main()
