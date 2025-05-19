import matplotlib.pyplot as plt
import numpy as np
import os
from constants import *

# Sample data
categories = ['Category A', 'Category B', 'Category C']
segment1 = [15, 25, 30]  # Bottom segment values
segment2 = [20, 15, 10]  # Middle segment values
segment3 = [10, 20, 15]  # Top segment values

# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the bottom segments
ax.bar(categories, segment1, label='Segment 1', color='skyblue')

# Plot middle segments on top of bottom ones
ax.bar(categories, segment2, bottom=segment1, label='Segment 2', color='salmon')

# Plot top segments on top of previous ones
ax.bar(categories, segment3, bottom=np.array(segment1)+np.array(segment2), 
       label='Segment 3', color='lightgreen')

# Add labels and title
ax.set_ylabel('Values')
ax.set_title('Stacked Bar Plot')
ax.legend()

plt.tight_layout()
os.makedirs(PLOTS, exist_ok=True)
plt.savefig(os.path.join(PLOTS, "test.png"))