#!/bin/bash
find . -iregex ".*\.\(cpp\|h\|cc\|hh\)" | etags  -
