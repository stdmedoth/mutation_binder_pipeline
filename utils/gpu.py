"""
GPU utility functions.
"""

import subprocess


def check_gpu_availability() -> bool:
    """
    Check if GPU is available (required for BioEmu).
    
    Returns:
        True if GPU is available, False otherwise
    """
    try:
        result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ GPU detected!")
            print(result.stdout)
            return True
        else:
            print("❌ No GPU detected.")
            return False
    except FileNotFoundError:
        print("❌ nvidia-smi not found. Make sure you're in a GPU-enabled environment.")
        return False


def get_gpu_info() -> dict:
    """
    Get detailed GPU information.
    
    Returns:
        Dictionary with GPU information or empty dict if no GPU
    """
    try:
        import torch
        if torch.cuda.is_available():
            return {
                'available': True,
                'device_count': torch.cuda.device_count(),
                'current_device': torch.cuda.current_device(),
                'device_name': torch.cuda.get_device_name(0),
                'memory_allocated': torch.cuda.memory_allocated(0),
                'memory_reserved': torch.cuda.memory_reserved(0),
            }
    except ImportError:
        pass
    
    return {'available': False}
