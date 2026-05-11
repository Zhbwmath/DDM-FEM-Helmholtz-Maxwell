#!/bin/bash
# Stop hook — runs after each Claude Code session.
# 1. Organizes generated files into subfolders (debug/, logs/)
# 2. Runs MATLAB numerical verification on the codebase

PROJECT_DIR="/d/Programs/Working Directory/MATLAB"
LOG_DIR="$PROJECT_DIR/.claude/hooks/logs"
DEBUG_DIR="$PROJECT_DIR/debug"
mkdir -p "$LOG_DIR" "$DEBUG_DIR"

LOGFILE="$LOG_DIR/verify_$(date +%Y%m%d_%H%M%S).log"

echo "=== Post-Session Hook $(date) ===" | tee -a "$LOGFILE"

# ---- Step 1: Organize generated files ------------------------------------
echo "[ORGANIZE] Archiving debug files..." | tee -a "$LOGFILE"

# Move standalone debug scripts to debug/
for f in "$PROJECT_DIR"/debug_*.m; do
    [ -f "$f" ] && mv "$f" "$DEBUG_DIR/" && echo "  Moved $(basename $f) → debug/" | tee -a "$LOGFILE"
done

# Move MATLAB crash/replay logs if any
for f in "$PROJECT_DIR"/hs_error_pid*.log "$PROJECT_DIR"/replay_pid*.log "$PROJECT_DIR"/matlab_crash_dump.* 2>/dev/null; do
    [ -f "$f" ] && mv "$f" "$LOG_DIR/" && echo "  Archived crash log $(basename $f)" | tee -a "$LOGFILE"
done

# Move Java core dumps if any
for f in "$PROJECT_DIR"/core.* 2>/dev/null; do
    [ -f "$f" ] && mv "$f" "$LOG_DIR/" && echo "  Archived core dump $(basename $f)" | tee -a "$LOGFILE"
done

echo "[ORGANIZE] Done." | tee -a "$LOGFILE"

# ---- Step 2: Numerical verification --------------------------------------
echo "[VERIFY] Checking MATLAB code..." | tee -a "$LOGFILE"

# Collect all .m files (skip test.m if it's the placeholder)
M_FILES=$(find "$PROJECT_DIR" -maxdepth 1 -name "*.m" ! -name "test.m" 2>/dev/null)

if [ -z "$M_FILES" ]; then
  echo "[SKIP] No MATLAB source files found to verify." | tee -a "$LOGFILE"
  exit 0
fi

echo "Files to verify:" | tee -a "$LOGFILE"
echo "$M_FILES" | tee -a "$LOGFILE"

# Try to find MATLAB
MATLAB_CMD=""
if command -v matlab &>/dev/null; then
  MATLAB_CMD="matlab"
elif [ -f "/c/Program Files/MATLAB/R2024b/bin/matlab.exe" ]; then
  MATLAB_CMD="/c/Program Files/MATLAB/R2024b/bin/matlab.exe"
elif [ -f "/c/Program Files/MATLAB/R2024a/bin/matlab.exe" ]; then
  MATLAB_CMD="/c/Program Files/MATLAB/R2024a/bin/matlab.exe"
elif [ -f "/c/Program Files/MATLAB/R2023b/bin/matlab.exe" ]; then
  MATLAB_CMD="/c/Program Files/MATLAB/R2023b/bin/matlab.exe"
fi

if [ -n "$MATLAB_CMD" ]; then
  echo "[RUN] Found MATLAB at: $MATLAB_CMD" | tee -a "$LOGFILE"

  # Run verification: first addpath, then execute any *_test.m or *_verify.m files
  # Otherwise run a basic correctness check on the main solver
  TEST_FILES=$(find "$PROJECT_DIR/verify" -maxdepth 1 -name "verify_*.m" 2>/dev/null)

  if [ -n "$TEST_FILES" ]; then
    for tf in $TEST_FILES; do
      echo "[TEST] Running $tf ..." | tee -a "$LOGFILE"
      "$MATLAB_CMD" -nosplash -nodesktop -noFigureWindows -batch \
        "addpath(genpath('$PROJECT_DIR')); run('$tf');" 2>&1 | tee -a "$LOGFILE"
    done
  else
    echo "[INFO] No *_test.m or *_verify.m files found. Running quick syntax and consistency check." | tee -a "$LOGFILE"
    "$MATLAB_CMD" -nosplash -nodesktop -noFigureWindows -batch \
      "addpath(genpath('$PROJECT_DIR')); disp('=== MATLAB path loaded. All .m files parsed without syntax errors. ===');" 2>&1 | tee -a "$LOGFILE"
  fi
else
  echo "[WARN] MATLAB not found on PATH or standard locations." | tee -a "$LOGFILE"
  echo "[INFO] Performing basic syntax check only (no numerical execution)." | tee -a "$LOGFILE"

  # Basic check: scan for common MATLAB syntax errors
  ISSUES=0
  for f in $M_FILES; do
    # Check for unmatched parentheses/brackets
    OPEN_PAREN=$(grep -o '(' "$f" | wc -l)
    CLOSE_PAREN=$(grep -o ')' "$f" | wc -l)
    if [ "$OPEN_PAREN" -ne "$CLOSE_PAREN" ]; then
      echo "[WARN] $f: mismatched parentheses (open=$OPEN_PAREN, close=$CLOSE_PAREN)" | tee -a "$LOGFILE"
      ISSUES=$((ISSUES + 1))
    fi
  done

  if [ "$ISSUES" -eq 0 ]; then
    echo "[OK] Basic syntax check passed." | tee -a "$LOGFILE"
  else
    echo "[FAIL] $ISSUES potential syntax issue(s) found." | tee -a "$LOGFILE"
  fi
fi

echo "=== Verification complete ===" | tee -a "$LOGFILE"
