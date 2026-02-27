#!/usr/bin/env python3
"""
检测特定细胞内的 Cy3 RNA 分子
使用 LoG 斑点检测识别小的球状 RNA 分子
统计细胞核内和细胞核外的 RNA 数量
"""

import argparse
import numpy as np
from skimage import io, feature, measure
from scipy import ndimage
import matplotlib.pyplot as plt
import json
import os
import czifile


def detect_rna_molecules_in_cell(rna_image, cell_mask, nucleus_mask, 
                                  min_sigma=0.5, max_sigma=3.0, 
                                  num_sigma=10, threshold=0.01):
    """
    在单个细胞内检测 RNA 分子
    
    参数:
        rna_image: 3D RNA 图像 (通道 0, Cy3)
        cell_mask: 细胞掩膜
        nucleus_mask: 细胞核掩膜
        min_sigma: LoG 最小 sigma
        max_sigma: LoG 最大 sigma
        num_sigma: sigma 数量
        threshold: 检测阈值
    
    返回:
        rna_coords: RNA 分子坐标列表 [(z, y, x), ...]
        rna_in_nucleus: 细胞核内 RNA 数量
        rna_in_cytoplasm: 细胞质内 RNA 数量
    """
    
    # 提取细胞内的 RNA 信号
    rna_in_cell = rna_image.copy()
    rna_in_cell[~cell_mask] = 0
    
    # 使用 LoG 斑点检测
    print(f"    运行 LoG 斑点检测...")
    print(f"      sigma 范围: {min_sigma}-{max_sigma}")
    print(f"      阈值: {threshold}")
    
    blobs = feature.blob_log(
        rna_in_cell,
        min_sigma=min_sigma,
        max_sigma=max_sigma,
        num_sigma=num_sigma,
        threshold=threshold,
        overlap=0.5
    )
    
    print(f"      检测到 {len(blobs)} 个候选 RNA 分子")
    
    if len(blobs) == 0:
        return [], 0, 0
    
    # 提取坐标 (z, y, x)
    rna_coords = blobs[:, :3].astype(int)
    
    # 过滤：确保在细胞内
    valid_rna = []
    for coord in rna_coords:
        z, y, x = coord
        if (0 <= z < cell_mask.shape[0] and 
            0 <= y < cell_mask.shape[1] and 
            0 <= x < cell_mask.shape[2]):
            if cell_mask[z, y, x]:
                valid_rna.append(coord)
    
    rna_coords = np.array(valid_rna)
    
    if len(rna_coords) == 0:
        return [], 0, 0
    
    # 统计细胞核内和细胞核外的 RNA
    rna_in_nucleus = 0
    rna_in_cytoplasm = 0
    
    for coord in rna_coords:
        z, y, x = coord
        if nucleus_mask[z, y, x]:
            rna_in_nucleus += 1
        else:
            rna_in_cytoplasm += 1
    
    return rna_coords, rna_in_nucleus, rna_in_cytoplasm


def create_rna_detection_report(rna_image, cell_labels, nuclei_mask, 
                                 cell_ids, rna_results, output_dir, z_slice=20):
    """
    创建 RNA 检测报告
    """
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'RNA Detection Report (Cells: {cell_ids})', fontsize=16, fontweight='bold')
    
    # 原始 RNA 图像 (Z=z_slice)
    ax = axes[0, 0]
    ax.imshow(rna_image[z_slice], cmap='hot', vmin=0, vmax=np.percentile(rna_image, 99))
    ax.set_title(f'RNA Signal (Cy3, Z={z_slice})')
    ax.axis('off')
    
    # 原始 RNA 图像 (MIP)
    ax = axes[0, 1]
    rna_mip = np.max(rna_image, axis=0)
    ax.imshow(rna_mip, cmap='hot', vmin=0, vmax=np.percentile(rna_mip, 99))
    ax.set_title('RNA Signal (Cy3, MIP)')
    ax.axis('off')
    
    # 细胞和细胞核 (Z=z_slice)
    ax = axes[0, 2]
    composite = np.zeros((*cell_labels.shape[1:], 3), dtype=np.uint8)
    for cell_id in cell_ids:
        cell_mask_2d = cell_labels[z_slice] == cell_id
        composite[cell_mask_2d] = [0, 255, 0]  # 绿色细胞
    nucleus_mask_2d = nuclei_mask[z_slice]
    composite[nucleus_mask_2d] = [0, 0, 255]  # 蓝色细胞核
    ax.imshow(composite)
    ax.set_title(f'Cells (green) + Nuclei (blue) (Z={z_slice})')
    ax.axis('off')
    
    # RNA 检测结果 (Z=z_slice)
    ax = axes[1, 0]
    ax.imshow(rna_image[z_slice], cmap='gray', vmin=0, vmax=np.percentile(rna_image, 99))
    for cell_id in cell_ids:
        if cell_id in rna_results:
            rna_coords = rna_results[cell_id]['coords']
            if len(rna_coords) > 0:
                # 只显示当前 Z 层的 RNA
                rna_in_slice = rna_coords[rna_coords[:, 0] == z_slice]
                if len(rna_in_slice) > 0:
                    ax.plot(rna_in_slice[:, 2], rna_in_slice[:, 1], 'r.', markersize=3)
    ax.set_title(f'RNA Molecules (red) (Z={z_slice})')
    ax.axis('off')
    
    # RNA 检测结果 (MIP)
    ax = axes[1, 1]
    ax.imshow(rna_mip, cmap='gray', vmin=0, vmax=np.percentile(rna_mip, 99))
    for cell_id in cell_ids:
        if cell_id in rna_results:
            rna_coords = rna_results[cell_id]['coords']
            if len(rna_coords) > 0:
                # MIP: 投影所有 Z 层
                ax.plot(rna_coords[:, 2], rna_coords[:, 1], 'r.', markersize=2, alpha=0.5)
    ax.set_title('RNA Molecules (red) (MIP)')
    ax.axis('off')
    
    # 统计信息
    ax = axes[1, 2]
    ax.axis('off')
    
    stats_text = "RNA Detection Statistics\n" + "="*40 + "\n\n"
    for cell_id in cell_ids:
        if cell_id in rna_results:
            result = rna_results[cell_id]
            total = result['total']
            if total > 0:
                stats_text += f"Cell {cell_id}:\n"
                stats_text += f"  Total RNA: {total}\n"
                stats_text += f"  In Nucleus: {result['in_nucleus']} ({result['in_nucleus']/total*100:.1f}%)\n"
                stats_text += f"  In Cytoplasm: {result['in_cytoplasm']} ({result['in_cytoplasm']/total*100:.1f}%)\n"
                stats_text += f"  Nucleus/Cytoplasm Ratio: {result['in_nucleus']/max(result['in_cytoplasm'], 1):.2f}\n\n"
            else:
                stats_text += f"Cell {cell_id}:\n"
                stats_text += f"  No RNA detected\n\n"
    
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    report_path = os.path.join(output_dir, 'rna_detection_report.png')
    plt.savefig(report_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  ✅ RNA 检测报告已保存: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description='检测特定细胞内的 Cy3 RNA 分子'
    )
    parser.add_argument('czi_file', help='CZI 文件路径')
    parser.add_argument('cell_labels', help='细胞分割标签文件')
    parser.add_argument('nuclei_mask', help='细胞核掩膜文件')
    parser.add_argument('cell_ids', nargs='+', type=int, help='要分析的细胞ID')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('--rna-channel', type=int, default=0, help='RNA 通道索引')
    parser.add_argument('--min-sigma', type=float, default=0.5, help='LoG 最小 sigma')
    parser.add_argument('--max-sigma', type=float, default=3.0, help='LoG 最大 sigma')
    parser.add_argument('--threshold', type=float, default=0.01, help='LoG 检测阈值')
    parser.add_argument('--z-slice', type=int, default=20, help='用于报告的 Z 切片')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("RNA 分子检测")
    print("=" * 60)
    
    print(f"\n读取 CZI 文件: {args.czi_file}")
    with czifile.CziFile(args.czi_file) as czi:
        image_data = czi.asarray()
    
    # 提取 RNA 通道
    print(f"  原始图像形状: {image_data.shape}")
    
    if image_data.ndim == 7:
        # (1, 1, C, Z, Y, X, 1)
        rna_image = image_data[0, 0, args.rna_channel, :, :, :, 0]
    elif image_data.ndim == 6:
        # (1, 1, C, Z, Y, X)
        rna_image = image_data[0, 0, args.rna_channel, :, :, :]
    else:
        raise ValueError(f"不支持的图像维度: {image_data.shape}")
    
    print(f"  RNA 图像形状: {rna_image.shape}")
    print(f"  RNA 图像范围: {rna_image.min()}-{rna_image.max()}")
    
    print(f"\n读取细胞分割标签: {args.cell_labels}")
    cell_labels = io.imread(args.cell_labels)
    print(f"  形状: {cell_labels.shape}")
    
    print(f"\n读取细胞核掩膜: {args.nuclei_mask}")
    nuclei_mask = io.imread(args.nuclei_mask)
    nuclei_mask = nuclei_mask > 0
    print(f"  形状: {nuclei_mask.shape}")
    
    print(f"\nLoG 斑点检测参数:")
    print(f"  min_sigma: {args.min_sigma}")
    print(f"  max_sigma: {args.max_sigma}")
    print(f"  threshold: {args.threshold}")
    
    # 检测每个细胞的 RNA
    rna_results = {}
    
    for cell_id in args.cell_ids:
        print(f"\n处理细胞 {cell_id}...")
        
        # 提取细胞和细胞核掩膜
        cell_mask = cell_labels == cell_id
        
        if not np.any(cell_mask):
            print(f"  ⚠️  细胞 {cell_id} 不存在")
            continue
        
        nucleus_mask_in_cell = nuclei_mask & cell_mask
        
        # 检测 RNA
        rna_coords, rna_in_nucleus, rna_in_cytoplasm = detect_rna_molecules_in_cell(
            rna_image, cell_mask, nucleus_mask_in_cell,
            min_sigma=args.min_sigma,
            max_sigma=args.max_sigma,
            threshold=args.threshold
        )
        
        total_rna = len(rna_coords)
        
        print(f"  检测到 {total_rna} 个 RNA 分子")
        if total_rna > 0:
            print(f"    细胞核内: {rna_in_nucleus} ({rna_in_nucleus/total_rna*100:.1f}%)")
            print(f"    细胞质内: {rna_in_cytoplasm} ({rna_in_cytoplasm/total_rna*100:.1f}%)")
        
        rna_results[cell_id] = {
            'coords': rna_coords,
            'total': total_rna,
            'in_nucleus': rna_in_nucleus,
            'in_cytoplasm': rna_in_cytoplasm
        }
        
        # 保存 RNA 坐标
        coords_file = os.path.join(args.output, f'cell_{cell_id}_rna_coords.npy')
        np.save(coords_file, rna_coords)
        print(f"  ✅ RNA 坐标已保存: {coords_file}")
    
    # 保存统计信息
    stats = {}
    for cell_id, result in rna_results.items():
        stats[str(cell_id)] = {
            'total': int(result['total']),
            'in_nucleus': int(result['in_nucleus']),
            'in_cytoplasm': int(result['in_cytoplasm']),
            'nucleus_cytoplasm_ratio': float(result['in_nucleus'] / max(result['in_cytoplasm'], 1))
        }
    
    stats_file = os.path.join(args.output, 'rna_statistics.json')
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"\n✅ 统计信息已保存: {stats_file}")
    
    # 创建检测报告
    print(f"\n创建 RNA 检测报告...")
    create_rna_detection_report(
        rna_image, cell_labels, nuclei_mask,
        args.cell_ids, rna_results, args.output,
        z_slice=args.z_slice
    )
    
    print("\n" + "=" * 60)
    print(f"✅ RNA 检测完成！共分析 {len(rna_results)} 个细胞")
    print("=" * 60)


if __name__ == '__main__':
    main()
