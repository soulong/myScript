#!/bin/bash
# Install dependencies for dna_barcode_new.py

echo "Installing Python dependencies for DNA barcode generator..."

# Try different methods to install pip
if ! python3 -m pip --version 2>/dev/null; then
    echo "pip not found, attempting to install..."
    
    # Try apt-get first
    if command -v apt-get &> /dev/null; then
        echo "Using apt-get to install python3-pip..."
        sudo apt-get update
        sudo apt-get install -y python3-pip python3-numpy python3-pandas python3-scipy python3-tqdm
    else
        echo "Please install pip first:"
        echo "  sudo apt-get install python3-pip"
        echo "Or use get-pip.py:"
        echo "  curl https://bootstrap.pypa.io/get-pip.py | python3"
        exit 1
    fi
fi

# Install remaining packages with pip
echo "Installing remaining packages..."
python3 -m pip install --user numpy pandas scipy tqdm

echo ""
echo "Dependencies installed successfully!"
echo "You can now run: python3 dna_barcode_new.py --length 16 --distance 4 --gc_min 40 --gc_max 60 --limit 100000 --cores 12"
