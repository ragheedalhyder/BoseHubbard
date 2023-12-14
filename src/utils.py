import os
import argparse
import yaml
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file", default="src/config.yml")
    args = parser.parse_args()
    return args


def create_output_dir():
    output_dir = os.path.join("res", datetime.now().strftime("%Y%m%d-%H%M%S"))
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def read_config(config_file):
    with open(config_file, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    return config
