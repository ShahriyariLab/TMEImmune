
FROM python:3.12-slim

WORKDIR /app

# Install required system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*


# Install required packages
RUN pip install --no-cache-dir TMEImmune numpy pandas scikit-learn cmapPy rnanorm gseapy inmoose lifelines matplotlib requests

# Define the default command
CMD ["python3"]
