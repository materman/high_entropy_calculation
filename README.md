![image](https://github.com/user-attachments/assets/777c0053-25fe-4acc-8b80-bc693cc02a9f)
![image](https://github.com/user-attachments/assets/0a373c3b-29aa-425f-860f-814192f93047)
执行命令  python calculate_sro.py POSCAR output.txt --cutoff 3.0
POSCAR：输入的结构文件路径。
output.txt：输出结果文件路径。
--cutoff：最近邻截断值（可选，默认3.0 Å）。
关键修改点
命令行参数解析
使用 argparse 模块解析输入文件、输出文件和截断值。
结果保存功能
将SRO参数矩阵和元素顺序写入指定输出文件。
提供默认截断值（3.0 Å），用户可通过 --cutoff 参数调整。
此脚本由deepseek生成。
