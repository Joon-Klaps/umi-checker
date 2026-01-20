#!/usr/bin/env bash
set -euo pipefail

REPO="Joon-Klaps/umi-checker"
BINARY_NAME="umi-checker"

OS=$(uname -s)
ARCH=$(uname -m)
TARGET_NAME="${BINARY_NAME}-${OS}-${ARCH}"

DEFAULT_INSTALL_DIR="$HOME/.local/bin"

LATEST_TAG=$(curl -fsSL "https://api.github.com/repos/${REPO}/releases/latest" \
  | grep '"tag_name":' | sed -E 's/.*"([^"]+)".*/\1/')

if [[ -z "${LATEST_TAG}" ]]; then
  echo "Could not determine latest release tag for ${REPO}."
  exit 1
fi

echo "Installing ${BINARY_NAME} ${LATEST_TAG} for ${OS}-${ARCH}"
echo

# --------------------------------------------------
# Choose install directory (interactive if possible)
# --------------------------------------------------
INSTALL_DIR="${DEFAULT_INSTALL_DIR}"

if [[ -t 0 ]]; then
  echo "Where do you want to install ${BINARY_NAME}?"
  echo "Press Enter to use default:"
  echo "  ${DEFAULT_INSTALL_DIR}"
  echo "Or enter a custom path (use '.' for current directory)"
  read -r -p "> " USER_INPUT

  if [[ -n "${USER_INPUT}" ]]; then
    if [[ "${USER_INPUT}" == "." ]]; then
      INSTALL_DIR="$(pwd)"
    else
      INSTALL_DIR="${USER_INPUT/#\~/$HOME}"
    fi
  fi
else
  echo "Non-interactive shell detected. Using default install dir: ${DEFAULT_INSTALL_DIR}"
fi

echo
echo "Installing to: ${INSTALL_DIR}"
echo

# --------------------------------------------------
# Download
# --------------------------------------------------
URL="https://github.com/${REPO}/releases/download/${LATEST_TAG}/${TARGET_NAME}"

HTTP_STATUS=$(curl -sSL -w "%{http_code}" -o "${BINARY_NAME}" "${URL}") || true
if [[ "${HTTP_STATUS}" != "200" ]]; then
  echo "Failed to download ${URL} (HTTP ${HTTP_STATUS})."
  echo "Available assets:"
  echo "  https://github.com/${REPO}/releases/tag/${LATEST_TAG}"
  exit 1
fi

chmod +x "${BINARY_NAME}"

mkdir -p "${INSTALL_DIR}"
mv "${BINARY_NAME}" "${INSTALL_DIR}/${BINARY_NAME}"

# --------------------------------------------------
# Final message
# --------------------------------------------------
cat <<EOF
--------------------------------------------------
Successfully installed ${BINARY_NAME} to:
  ${INSTALL_DIR}/${BINARY_NAME}

EOF

if [[ ":$PATH:" != *":${INSTALL_DIR}:"* ]]; then
  cat <<EOF
Note: ${INSTALL_DIR} is not in your PATH.
You may want to add it:

  export PATH="${INSTALL_DIR}:\$PATH"

EOF
fi

echo "Run: ${BINARY_NAME} --version"
