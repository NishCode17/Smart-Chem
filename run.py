import torch

if torch.cuda.is_available():
    print("🚀 good news: Found NVIDIA GPU!")
    print(f"GPU Name: {torch.cuda.get_device_name(0)}")
elif torch.backends.mps.is_available():
    print("🍎 Good news: Found Apple Silicon (M1/M2/M3)!")
else:
    print("🐌 Warning: No GPU found. Training will run on CPU (Slow but feasible).")