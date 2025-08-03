#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
# import dash
# # from dash.exceptions import PreventUpdate
# from dash import dcc, html
# from dash.dependencies import Input, Output
# from flask import Flask
# import plotly.express as px
# import pandas
# import pickle

# import os
# os.chdir('/media/hao/Data1/single_cell_dataset/hpa-single-cell-image-classification')
import sys
sys.path.append('/home/hao/Documents/GitHub/myScript/python_functions')
from dl_helper import vis_feature


# Run Dash App
if __name__ == "__main__":
    
    # Argument parsing
    parser = argparse.ArgumentParser(description="Interactive UMAP Visualization")
    parser.add_argument("--data_path", type=str, required=True, help="Path to the pickle file containing embeddings.")
    args = parser.parse_args()
    
    vis_feature(args.data_path)



