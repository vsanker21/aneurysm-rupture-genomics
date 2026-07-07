"""
Install All Required Packages for Advanced ML Analyses
"""
import subprocess
import sys

print("=" * 80)
print("INSTALLING REQUIRED PACKAGES")
print("=" * 80)

packages = [
    'tensorflow>=2.10.0',
    'scikit-learn>=1.0.0',
    'pandas>=1.3.0',
    'numpy>=1.21.0',
    'matplotlib>=3.5.0',
    'seaborn>=0.11.0',
    'requests>=2.28.0'
]

print("\nInstalling packages...")
for package in packages:
    print(f"\nInstalling {package}...")
    try:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', package, '--quiet'])
        print(f"  [OK] {package} installed/updated")
    except subprocess.CalledProcessError as e:
        print(f"  [ERROR] Failed to install {package}: {e}")
        # Try without version requirement
        try:
            pkg_name = package.split('>=')[0]
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', pkg_name, '--quiet'])
            print(f"  [OK] {pkg_name} installed (without version check)")
        except:
            print(f"  [WARNING] Could not install {pkg_name}")

print("\n" + "=" * 80)
print("VERIFYING INSTALLATIONS")
print("=" * 80)

imported = []
failed = []

for package in ['tensorflow', 'sklearn', 'pandas', 'numpy', 'matplotlib', 'seaborn', 'requests']:
    try:
        __import__(package.replace('-', '_'))
        imported.append(package)
        print(f"[OK] {package} imported successfully")
    except ImportError:
        failed.append(package)
        print(f"[ERROR] {package} import failed")

if failed:
    print(f"\n[WARNING] Failed to import: {', '.join(failed)}")
    print("Please install manually: pip install " + " ".join(failed))
else:
    print("\n[OK] All required packages installed and verified!")

print("=" * 80)

