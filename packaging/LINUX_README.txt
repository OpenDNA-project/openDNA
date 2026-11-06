openDNA for Linux
=================

Run ./openDNA in this folder. Keep the folder together: the program needs the _internal folder next to it.
Free software (AGPL-3.0). Download it only from https://github.com/OpenDNA-project/openDNA/releases . If you paid for it, you were scammed.

Built on Ubuntu 22.04 (x86-64): it runs on that and on newer distributions (glibc 2.35 or later). It needs a desktop (X11 or
Wayland) with the usual graphics libraries (libGL / libEGL, fontconfig); Qt and the X11 helper libraries it uses are inside.
If it does not start, run it from a terminal to see why; on Debian or Ubuntu
    sudo apt install libegl1 libgl1 libxkbcommon-x11-0 libxcb-cursor0 libfontconfig1
installs what a minimal system may lack.

Read mapping (Tools > Install Mapping Tools) runs bwa, samtools and minimap2 directly (no WSL as on Windows).

Using another Qt (the LGPL-3.0 lets you)
----------------------------------------
openDNA uses Qt through PySide6. In this folder Qt is made of plain files, so you can swap them: get the PySide6-Essentials and
shiboken6 wheels of the version you want for Linux x86-64, Python 3.13 (see THIRD_PARTY_NOTICES.md), close openDNA, rename the
folders _internal/PySide6 and _internal/shiboken6, unzip each wheel and copy its PySide6 and shiboken6 folders into _internal.
Help > About shows the Qt version in use.

This version does not update itself: download the new one from the releases page when there is one.

Licences and the source of everything in here: THIRD_PARTY_NOTICES.md and LICENSE (next to this file). Source code: https://github.com/OpenDNA-project/openDNA
