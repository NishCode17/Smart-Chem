
import os
import subprocess
import time
import sys

def force_kill_port_8000():
    print("🔍 Scanning for processes on port 8000...")
    
    # Get all PIDs on port 8000
    cmd = "netstat -ano | findstr :8000"
    try:
        output = subprocess.check_output(cmd, shell=True).decode()
        lines = output.strip().split('\n')
        
        pids_to_kill = set()
        for line in lines:
            if "LISTENING" in line:
                parts = line.split()
                pid = parts[-1]
                if pid != "0":
                    pids_to_kill.add(pid)
                    
        if not pids_to_kill:
            print("✅ Port 8000 is free (No LISTENING process found).")
            return

        print(f"⚠️ Found processes holding stuck port: {pids_to_kill}")
        
        for pid in pids_to_kill:
            print(f"💥 Killing PID {pid}...")
            os.system(f"taskkill /F /PID {pid}")
            
        print("✅ Cleanup complete. Port 8000 should be free.")
        
    except subprocess.CalledProcessError:
        print("✅ Port 8000 seems free (netstat found nothing).")
    except Exception as e:
        print(f"❌ Error during cleanup: {e}")

if __name__ == "__main__":
    force_kill_port_8000()
