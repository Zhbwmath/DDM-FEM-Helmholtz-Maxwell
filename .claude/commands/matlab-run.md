---
name: matlab-run
description: Run MATLAB code or scripts silently (no console, figures allowed). Usage: /matlab-run <script.m> or /matlab-run "<MATLAB code>"
---

# MATLAB Silent Runner

Run MATLAB code without showing the desktop/console. Figure windows are allowed to display.

## Usage
- `/matlab-run verify/verify_all.m` — run a script
- `/matlab-run "assembleNed2CurlCurl2D; disp('done');"` — run inline code
- `/matlab-run "mesh = squaremesh([0 1 0 1], 0.1); spy(assembleCurlCurl2D(mesh.node, mesh.elem));"` — run and show a figure

## Execution

Determine the MATLAB executable:
- Check `C:\Program Files\MATLAB\R2024b\bin\matlab.exe` first
- Then `C:\Program Files\MATLAB\R2024a\bin\matlab.exe`
- Fall back to `matlab` on PATH

Run with flags: `-nosplash -nodesktop -batch "<code>"`

The `-batch` flag:
- Starts MATLAB without the desktop (no console clutter)
- Executes the given statement and exits
- **Allows figure windows** to appear (unlike `-noFigureWindows`)
- Returns non-zero exit code on error

Always prepend `addpath(genpath('.'));` so all project functions are on the path.

## Examples

```bash
# Run a verification script
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); run('verify/verify_ned1_2D.m');"

# Run inline code with figure
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); mesh=cubemesh([0 1 0 1 0 1],0.2); A=assembleCurlCurl3D(mesh.node,mesh.elem); spy(A);"

# Quick computation
matlab -nosplash -nodesktop -batch "addpath(genpath('.')); mesh=squaremesh([0 1 0 1],0.1); A=assembleCurlCurl2D(mesh.node,mesh.elem); [V,D]=eigs(A,6,'smallestabs'); disp(diag(D));"
```
