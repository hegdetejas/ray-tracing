import numpy as np
import matplotlib.pyplot as plt
import ast


def extract_data(filename):
    f = open(filename, "r")

    file_str = f.read()
    data = ast.literal_eval(file_str)
    NT = data[0][0]
    grid = data[1]
    np_grid = np.array(grid)

    return (NT, np_grid)


def visualize(filename, NT, np_grid):
    plt.figure(figsize=(10, 8))
    plt.imshow(np_grid, cmap=plt.cm.get_cmap("binary_r"))
    plt.savefig(filename)
    plt.clf()


if __name__ == "__main__":
    data = extract_data("../out/sphere.txt")
    visualize("../out/sphere.png", data[0], data[1])
