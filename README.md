# Transportation Optimization (VAM & MODI)

This repository contains a Python implementation of the Transportation Problem, a classic Operations Research optimization challenge. The script calculates an Initial Basic Feasible Solution using Vogel's Approximation Method (VAM) and optimizes it using the Modified Distribution (MODI) method.

## Features
* Validates supply and demand to ensure a balanced problem.
* Automates VAM penalty calculations and allocations.
* Performs MODI optimality testing and generates improvement loops.
* Outputs step-by-step terminal visualizations of the cost matrices.

## How to Run
1. Ensure you have Python 3.x installed.
2. Clone this repository.
3. Run the script from your terminal:
   `python main.py`

## Input Format
The program expects the following inputs via the terminal:
* Number of sources and destinations
* The transportation cost matrix
* Supply values for each source
* Demand values for each destination

## Example Output
<img width="848" height="510" alt="Screenshot 2026-09-13 191308" src="https://github.com/user-attachments/assets/6072e8bd-cb16-47a9-8e9d-efbedd39f5cd" />
<img width="843" height="519" alt="Screenshot 2026-09-13 191340" src="https://github.com/user-attachments/assets/7279535c-f847-410c-986c-4836b18e6dc1" />

