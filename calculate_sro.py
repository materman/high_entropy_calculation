mport numpy as np
from scipy.spatial import KDTree
import argparse

def parse_poscar(poscar_path):
    """
    解析POSCAR文件，返回晶格矩阵、元素列表、原子坐标（笛卡尔坐标）
    """
    with open(poscar_path, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    # 提取晶格缩放因子和晶格向量
    scale_factor = float(lines[1])
    lattice = np.array([list(map(float, lines[i].split())) for i in range(2,5)]) * scale_factor

    # 提取元素种类和原子数量
    elements = lines[5].split()
    atom_counts = list(map(int, lines[6].split()))

    # 确定坐标类型（Direct或Cartesian）
    coord_type = lines[7].lower()

    # 提取原子坐标并转换为笛卡尔坐标
    coords = []
    species = []
    idx = 8
    for elem, count in zip(elements, atom_counts):
        for _ in range(count):
            parts = lines[idx].split()
            frac_coord = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
            if coord_type == 'direct':
                # 分数坐标转笛卡尔坐标: [x, y, z] @ lattice_matrix
                cart_coord = np.dot(frac_coord, lattice)
            else:
                cart_coord = frac_coord * scale_factor
            coords.append(cart_coord)
            species.append(elem)
            idx += 1

    return lattice, species, np.array(coords)

def calculate_sro(poscar_path, cutoff):
    """
    计算Warren-Cowley SRO参数（仅依赖scipy和numpy）
    """
    # 解析POSCAR
    lattice, species, coords = parse_poscar(poscar_path)
    elements = list(sorted(set(species)))
    n_elements = len(elements)
    X = {elem: species.count(elem)/len(species) for elem in elements}

    # 构建KDTree（考虑周期性边界）
    # 计算晶格各方向长度（用于周期性镜像）
    a, b, c = np.linalg.norm(lattice, axis=1)
    kdtree = KDTree(coords, boxsize=[a, b, c])  # 关键：设置周期性边界

    # 统计原子对
    pair_counts = np.zeros((n_elements, n_elements), dtype=int)
    for i in range(len(coords)):
        # 查找i原子的最近邻（包含自身，需排除）
        neighbors = kdtree.query_ball_point(coords[i], cutoff)
        neighbors = [j for j in neighbors if j != i]  # 移除自身
        # 统计元素对
        elem_i = elements.index(species[i])
        for j in neighbors:
            elem_j = elements.index(species[j])
            pair_counts[elem_i, elem_j] += 1

    # 对称化矩阵（避免重复计数）
    pair_counts = (pair_counts + pair_counts.T) // 2

    # 计算SRO参数
    total_pairs = pair_counts.sum()
    sro_matrix = np.zeros((n_elements, n_elements))
    for i in range(n_elements):
        for j in range(n_elements):
            X_i = X[elements[i]]
            X_j = X[elements[j]]
            if X_i * X_j == 0:
                sro = 0.0
            else:
                P_ij = pair_counts[i, j] / total_pairs
                sro = 1 - (P_ij / (X_i * X_j))
            sro_matrix[i, j] = sro

    return sro_matrix, elements

def save_sro_results(output_path, sro_matrix, elements):
    """
    将SRO参数矩阵保存到文件
    """
    with open(output_path, 'w') as f:
        # 写入元素顺序
        f.write("元素顺序: " + "  ".join(elements) + "\n")
        # 写入SRO矩阵
        f.write("Warren-Cowley SRO参数矩阵:\n")
        for i in range(len(elements)):
            row = [f"{sro_matrix[i,j]:.3f}" for j in range(len(elements))]
            f.write("  " + "  ".join(row) + "\n")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="计算高熵合金的Warren-Cowley SRO参数")
    parser.add_argument("input_file", help="输入POSCAR文件路径")
    parser.add_argument("output_file", help="输出结果文件路径")
    parser.add_argument("--cutoff", type=float, default=3.0, help="最近邻截断值（默认：3.0 Å）")
    args = parser.parse_args()

    # 计算SRO参数
    sro_matrix, elements = calculate_sro(args.input_file, args.cutoff)

    # 保存结果
    save_sro_results(args.output_file, sro_matrix, elements)
    print(f"结果已保存到 {args.output_file}")

if __name__ == "__main__":
    main()
