#!/bin/bash

OS=$(uname -s)
# Define the .zshrc (.bashrc) file path
ZSHRC_FILE="$HOME/.zshrc"

if [ "$OS" == "Darwin" ]; then
  echo "-> Operating System Detected: macOS ..."
  ZSHRC_FILE="${HOME}/.zshrc" # In macOS default shell is zsh hence zshrc
  # Get macOS version
  MACOS_VERSION=$(sw_vers -productVersion)
  echo "-> macOS Version: ${MACOS_VERSION} ..."
fi

if [ "$OS" == "Linux" ]; then
  echo "-> Operating System Detected: Linux ..."
  if [[ "$SHELL" == *"bash"* ]]; then
    ZSHRC_FILE="${HOME}/.bashrc" # In Linux default shell is
  elif [[ "$SHELL" == *"zsh"* ]]; then
    ZSHRC_FILE="${HOME}/.bashrc" # In Linux zsh shell is being used
  else
    echo "-> Unknown shell !!."
    exit 1
  fi
  # Get kernel version
  KERNEL_VERSION=$(uname -r)
  echo "-> Kernel Version: ${KERNEL_VERSION} ..."

  # Get distribution details from /etc/os-release
  if [ -f /etc/os-release ]; then
    . /etc/os-release
    DISTRIBUTION_NAME=$NAME
    DISTRIBUTION_VERSION=$VERSION_ID
    echo "-> Distribution: ${DISTRIBUTION_NAME} ..."
    echo "-> Distribution Version: ${DISTRIBUTION_VERSION} ..."
  else
    echo "-> Distribution details not found ..."
  fi
fi

#sleep for 2 seconds
sleep 1

# Define the repository URL
REPO_URL="https://github.com/luketpickering/ROOTEGPythia6.git"

# Define the target directory name
DIR_NAME="ROOTEGPythia6"

# Define the full path for the target directory
TARGET_DIR="$HOME/$DIR_NAME"

# Print the default path and ask for confirmation or a new path
echo "-> ROOTEGPythia6 will be cloned into ${DIR_NAME} directory ..."
echo "->* credit: Luke Pickering, for more details please visit https://github.com/luketpickering/ROOTEGPythia6.git ..."
sleep 2
echo "-> ${DIR_NAME} will be created in ${HOME}/ directory."
read -p "-> Do you want to proceed with the default path (y) or provide a new path? [y/n]: " choice

if [[ "$choice" =~ ^[Nn]$ ]]; then
  # Prompt the user for a new path
  read -p "-> Please enter the new absolute path for ${DIR_NAME}: " new_path

  # Check if the entered path is valid
  if [ -d "$new_path" ]; then
    TARGET_DIR="${new_path}/${DIR_NAME}"
    echo "-> The new path for ${DIR_NAME} will be: ${TARGET_DIR} ..."
  else
    echo "-> Error: The provided path '${new_path}' is not a valid directory ..."
    exit 1
  fi
else
  echo "-> Proceeding with the default path: ${TARGET_DIR} ..."
fi

# Check if the directory already exists
if [ -d "$TARGET_DIR" ]; then
  echo "-> Directory $TARGET_DIR already exists... Entering $TARGET_DIR ..."
  sleep 2
  cd "$TARGET_DIR" || exit 1 # '|| exit 1' exits the script if the cd command fails
else
  # Create the directory
  echo "-> Creating directory $TARGET_DIR..."
  mkdir -p "$TARGET_DIR"

  # Clone the repository into the new directory
  git clone "$REPO_URL" "$TARGET_DIR"

  # Check if the clone was successful
  if [ $? -eq 0 ]; then
    echo "-> Repository cloned successfully. Entering $TARGET_DIR ..."
    cd "$TARGET_DIR" || exit 1 # '|| exit 1' exits the script if the cd command fails
  else
    echo "** Error: Failed to clone the repository ..."
    exit 1
  fi
fi

# Check if the directory already exists
if [ -d "build" ]; then
  echo "-> build/ directory already exists. Emptying build/ directory..."
  rm -r build
  mkdir -p build && cd build || exit 1 # '&&' ensures the second command only runs if the first succeeds
else
  # Make the 'build' directory and jump into it
  echo "-> Making and entering 'build/' directory..."
  mkdir -p build && cd build || exit 1 # '&&' ensures the second command only runs if the first succeeds
fi

# Run CMake and make install
echo "-> Running CMake and make install..."
cmake .. -DROOTEGPythia6_Pythia6_BUILTIN=ON && make install -j 4

# Check if the commands were successful
if [ $? -eq 0 ]; then
  echo "-> Build and installation completed successfully ..."
else
  echo "** Error: Build or installation failed ..."
  exit 1
fi

echo "-> ROOTEGPythia6 installed successfully ..."

# Define the variable name and its value
ENVM_VAR_NAME="ROOTEGPythia6_ROOT"
ENVM_VAR_VALUE="${TARGET_DIR}/build/${OS}"

echo "-> Checking whether environment variable ${ENVM_VAR_NAME} exists ..."
sleep 1

# Check if the variable is already in .zshrc
if ! grep -q "export $ENVM_VAR_NAME=" "$ZSHRC_FILE"; then
  # Append the export command to the .zshrc file
  echo "export ${ENVM_VAR_NAME}=\"${ENVM_VAR_VALUE}\"" >>"${ZSHRC_FILE}"
  echo "-> Added ${ENVM_VAR_NAME} to ${ZSHRC_FILE} ..."
else
  echo "** Environment variable ${ENVM_VAR_NAME} already exists in ${ZSHRC_FILE} ..."
fi

sleep 1
echo "-> Please source your ${ZSHRC_FILE} file to load changes ..."
echo "-> After that you can now prceed with NuWro installation."
