# VMD stdout redirector
# 
# Author: Robin Betz
# 
# Copyright (C) 2015 Robin Betz
# 
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330
# Boston, MA 02111-1307, USA.

import os, sys

class VmdSilencer:
    """
    Toggles whether or not C extensions can write to stdout. Since
    VMD is the only C extension that does this, silences this extra
    output during the dabbling process. This is done in dabble.py not
    in the Dabble API because the user probably wants more info / it 
    can't hurt them if they're smart enough to use the API.

    Most of this is from:
    http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/
    which is licensed under the MIT license.

    Attributes:
        output (file): Where to put the VMD ouptut
    """

    def __init__(self, output=os.devnull):
        self.outfile = output
        self.mode = 'w'
        
    #==========================================================================

    def __enter__(self):
        self.sys = sys
        # save previous stdout/stderr
        self.saved_stream = sys.__stdout__
        self.fd = self.saved_stream.fileno()
        self.saved_fd = os.dup(sys.stdout.fileno())
        sys.stdout.flush() # flush any pending output 

        # Check if we're not actually redirecting
        if self.outfile == sys.stdout:
            return

        # open surrogate files
        null_fd = open(self.outfile, self.mode)
        os.dup2(null_fd.fileno(), self.fd)

        self.null_stream = open(self.outfile, self.mode, 0)
        self.null_fd = self.null_stream.fileno()

        # overwrite file objects and low-level file descriptors
        os.dup2(self.null_fd, self.fd)

        sys.stdout = os.fdopen(self.saved_fd, 'w')
        sys.stdout.flush()

    #==========================================================================

    def __exit__(self, *args):
        sys = self.sys
        # flush any pending output
        sys.__stdout__.flush()
        # restore original streams and file descriptors
        os.dup2(self.saved_fd, self.fd)
        sys.stdout = self.saved_stream
        # clean up
        try: self.null_stream.close()
        except: pass
        os.close(self.saved_fd)
        return False

    #==========================================================================

