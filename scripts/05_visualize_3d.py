#!/usr/bin/env python3
"""
为单个细胞创建包含 RNA 分子的完整 3D 可视化（修正版）
使用平滑后的细胞轮廓过滤 RNA 分子，确保只显示轮廓内的 RNA
包括：细胞轮廓（平滑表面）+ 细胞核（平滑球形）+ RNA 分子（小球）
"""

import argparse
import numpy as np
from skimage import io, measure
from scipy import ndimage
import plotly.graph_objects as go
import json
import os


def smooth_mask(mask, sigma=1.0):
    """对掩膜进行高斯平滑"""
    mask_float = mask.astype(float)
    smoothed = ndimage.gaussian_filter(mask_float, sigma=sigma)
    return smoothed


def create_smoothed_mask(mask, sigma=1.0, threshold=0.5):
    """
    创建平滑后的二值掩膜
    
    参数:
        mask: 原始二值掩膜
        sigma: 高斯平滑参数
        threshold: 平滑后的阈值
    
    返回:
        smoothed_binary_mask: 平滑后的二值掩膜
    """
    smoothed = smooth_mask(mask, sigma=sigma)
    smoothed_binary = smoothed >= threshold
    return smoothed_binary


def create_mesh_from_mask(mask, level=0.5, smooth_sigma=1.0):
    """从掩膜创建 3D 网格，使用 Marching Cubes 算法"""
    smoothed_mask = smooth_mask(mask, sigma=smooth_sigma)
    
    try:
        verts, faces, normals, values = measure.marching_cubes(
            smoothed_mask, 
            level=level,
            spacing=(1.0, 1.0, 1.0)
        )
        return verts, faces
    except Exception as e:
        print(f"    ⚠️  Marching Cubes 失败: {e}")
        return None, None


def filter_rna_by_smoothed_mask(rna_coords, smoothed_cell_mask, smoothed_nucleus_mask):
    """
    根据平滑后的掩膜过滤 RNA 分子
    
    参数:
        rna_coords: RNA 分子坐标 (N, 3) 数组 [(z, y, x), ...]
        smoothed_cell_mask: 平滑后的细胞掩膜
        smoothed_nucleus_mask: 平滑后的细胞核掩膜
    
    返回:
        filtered_rna_coords: 过滤后的 RNA 坐标
        rna_in_nucleus: 细胞核内 RNA 数量
        rna_in_cytoplasm: 细胞质内 RNA 数量
    """
    
    if len(rna_coords) == 0:
        return np.array([]), 0, 0
    
    # 过滤：只保留在平滑后细胞掩膜内的 RNA
    valid_rna = []
    rna_in_nucleus = 0
    rna_in_cytoplasm = 0
    
    for coord in rna_coords:
        z, y, x = coord.astype(int)
        
        # 检查坐标是否在图像范围内
        if (0 <= z < smoothed_cell_mask.shape[0] and 
            0 <= y < smoothed_cell_mask.shape[1] and 
            0 <= x < smoothed_cell_mask.shape[2]):
            
            # 检查是否在平滑后的细胞内
            if smoothed_cell_mask[z, y, x]:
                valid_rna.append(coord)
                
                # 统计细胞核内/外
                if smoothed_nucleus_mask[z, y, x]:
                    rna_in_nucleus += 1
                else:
                    rna_in_cytoplasm += 1
    
    filtered_rna_coords = np.array(valid_rna) if valid_rna else np.array([])
    
    return filtered_rna_coords, rna_in_nucleus, rna_in_cytoplasm


def create_3d_visualization_with_rna(cell_labels, nuclei_mask, rna_coords, cell_id, 
                                      output_path, cell_smooth_sigma=1.5, 
                                      nucleus_smooth_sigma=1.0, rna_sphere_size=1.0):
    """
    为单个细胞创建包含 RNA 分子的完整 3D 可视化（修正版）
    
    参数:
        cell_labels: 细胞分割标签
        nuclei_mask: 细胞核掩膜
        rna_coords: RNA 分子坐标 (N, 3) 数组 [(z, y, x), ...]
        cell_id: 细胞 ID
        output_path: 输出文件路径
        cell_smooth_sigma: 细胞轮廓平滑参数
        nucleus_smooth_sigma: 细胞核平滑参数
        rna_sphere_size: RNA 小球大小
    """
    
    print(f"\n为细胞 {cell_id} 创建包含 RNA 的 3D 可视化（修正版）...")
    
    # 提取细胞掩膜
    cell_mask = cell_labels == cell_id
    
    if not np.any(cell_mask):
        print(f"  ⚠️  细胞 {cell_id} 不存在")
        return
    
    # 提取细胞内的细胞核
    nucleus_in_cell = nuclei_mask & cell_mask
    has_nucleus = np.any(nucleus_in_cell)
    
    print(f"  原始细胞体积: {np.sum(cell_mask)} px³")
    if has_nucleus:
        print(f"  原始细胞核体积: {np.sum(nucleus_in_cell)} px³")
    print(f"  原始 RNA 分子数量: {len(rna_coords)}")
    
    # 创建平滑后的掩膜
    print(f"  创建平滑后的掩膜...")
    smoothed_cell_mask = create_smoothed_mask(cell_mask, sigma=cell_smooth_sigma, threshold=0.5)
    smoothed_nucleus_mask = create_smoothed_mask(nucleus_in_cell, sigma=nucleus_smooth_sigma, threshold=0.5)
    
    print(f"  平滑后细胞体积: {np.sum(smoothed_cell_mask)} px³")
    if has_nucleus:
        print(f"  平滑后细胞核体积: {np.sum(smoothed_nucleus_mask)} px³")
    
    # 根据平滑后的掩膜过滤 RNA
    print(f"  根据平滑后的掩膜过滤 RNA...")
    filtered_rna_coords, rna_in_nucleus, rna_in_cytoplasm = filter_rna_by_smoothed_mask(
        rna_coords, smoothed_cell_mask, smoothed_nucleus_mask
    )
    
    total_rna = len(filtered_rna_coords)
    print(f"  过滤后 RNA 分子数量: {total_rna}")
    if total_rna > 0:
        print(f"    细胞核内: {rna_in_nucleus} ({rna_in_nucleus/total_rna*100:.1f}%)")
        print(f"    细胞质内: {rna_in_cytoplasm} ({rna_in_cytoplasm/total_rna*100:.1f}%)")
    
    # 创建 3D 图形
    fig = go.Figure()
    
    # 1. 生成细胞轮廓网格（绿色，半透明）
    print(f"  生成细胞轮廓网格...")
    cell_verts, cell_faces = create_mesh_from_mask(cell_mask, level=0.5, smooth_sigma=cell_smooth_sigma)
    
    if cell_verts is not None and cell_faces is not None:
        fig.add_trace(go.Mesh3d(
            x=cell_verts[:, 2],
            y=cell_verts[:, 1],
            z=cell_verts[:, 0],
            i=cell_faces[:, 0],
            j=cell_faces[:, 1],
            k=cell_faces[:, 2],
            color='lightgreen',
            opacity=0.2,
            name='Cell Boundary (CFW)',
            flatshading=False,
            lighting=dict(
                ambient=0.5,
                diffuse=0.8,
                specular=0.2,
                roughness=0.5
            ),
            lightposition=dict(x=100, y=100, z=100)
        ))
        print(f"    ✓ 细胞轮廓: {len(cell_verts)} 顶点, {len(cell_faces)} 面")
    
    # 2. 生成细胞核网格（蓝色，半透明）
    if has_nucleus:
        print(f"  生成细胞核网格...")
        nucleus_verts, nucleus_faces = create_mesh_from_mask(nucleus_in_cell, level=0.5, smooth_sigma=nucleus_smooth_sigma)
        
        if nucleus_verts is not None and nucleus_faces is not None:
            fig.add_trace(go.Mesh3d(
                x=nucleus_verts[:, 2],
                y=nucleus_verts[:, 1],
                z=nucleus_verts[:, 0],
                i=nucleus_faces[:, 0],
                j=nucleus_faces[:, 1],
                k=nucleus_faces[:, 2],
                color='blue',
                opacity=0.6,
                name='Nucleus (DAPI)',
                flatshading=False,
                lighting=dict(
                    ambient=0.4,
                    diffuse=0.9,
                    specular=0.5,
                    roughness=0.3
                ),
                lightposition=dict(x=100, y=100, z=100)
            ))
            print(f"    ✓ 细胞核: {len(nucleus_verts)} 顶点, {len(nucleus_faces)} 面")
    
    # 3. 添加 RNA 分子（红色小球）
    if len(filtered_rna_coords) > 0:
        print(f"  添加 RNA 分子...")
        
        # 使用 Scatter3d 绘制 RNA 分子为小球
        fig.add_trace(go.Scatter3d(
            x=filtered_rna_coords[:, 2],
            y=filtered_rna_coords[:, 1],
            z=filtered_rna_coords[:, 0],
            mode='markers',
            marker=dict(
                size=rna_sphere_size,
                color='red',
                opacity=0.8,
                symbol='circle'
            ),
            name=f'RNA Molecules (Cy3, N={len(filtered_rna_coords)})'
        ))
        print(f"    ✓ RNA 分子: {len(filtered_rna_coords)} 个")
    
    # 设置布局
    fig.update_layout(
        title=f'3D Visualization of Cell {cell_id} with RNA Molecules (Corrected)',
        scene=dict(
            xaxis_title='X (px)',
            yaxis_title='Y (px)',
            zaxis_title='Z (px)',
            aspectmode='data',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5)
            )
        ),
        width=1200,
        height=900,
        showlegend=True
    )
    
    # 保存为 HTML
    fig.write_html(output_path)
    print(f"  ✅ 3D 可视化已保存: {output_path}")
    
    # 返回统计信息
    return {
        'cell_id': cell_id,
        'original_rna_count': len(rna_coords),
        'filtered_rna_count': total_rna,
        'rna_in_nucleus': rna_in_nucleus,
        'rna_in_cytoplasm': rna_in_cytoplasm
    }


def main():
    parser = argparse.ArgumentParser(
        description='为单个细胞创建包含 RNA 分子的完整 3D 可视化（修正版）'
    )
    parser.add_argument('cell_labels', help='细胞分割标签文件')
    parser.add_argument('nuclei_mask', help='细胞核掩膜文件')
    parser.add_argument('rna_dir', help='RNA 坐标文件目录')
    parser.add_argument('cell_ids', nargs='+', type=int, help='要可视化的细胞ID')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('--cell-smooth', type=float, default=1.5, 
                        help='细胞轮廓平滑参数 (sigma)')
    parser.add_argument('--nucleus-smooth', type=float, default=1.0, 
                        help='细胞核平滑参数 (sigma)')
    parser.add_argument('--rna-size', type=float, default=2.0, 
                        help='RNA 小球大小')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("单个细胞完整 3D 可视化（含 RNA 分子，修正版）")
    print("=" * 60)
    
    print(f"\n读取细胞分割标签: {args.cell_labels}")
    cell_labels = io.imread(args.cell_labels)
    print(f"  形状: {cell_labels.shape}")
    
    print(f"\n读取细胞核掩膜: {args.nuclei_mask}")
    nuclei_mask = io.imread(args.nuclei_mask)
    nuclei_mask = nuclei_mask > 0
    print(f"  形状: {nuclei_mask.shape}")
    
    print(f"\n可视化参数:")
    print(f"  细胞轮廓 sigma: {args.cell_smooth}")
    print(f"  细胞核 sigma: {args.nucleus_smooth}")
    print(f"  RNA 小球大小: {args.rna_size}")
    
    # 为每个细胞创建 3D 可视化
    all_stats = {}
    
    for cell_id in args.cell_ids:
        # 读取 RNA 坐标
        rna_coords_file = os.path.join(args.rna_dir, f'cell_{cell_id}_rna_coords.npy')
        
        if not os.path.exists(rna_coords_file):
            print(f"\n⚠️  细胞 {cell_id} 的 RNA 坐标文件不存在: {rna_coords_file}")
            continue
        
        rna_coords = np.load(rna_coords_file)
        
        output_path = os.path.join(args.output, f'cell_{cell_id}_with_rna_3d_corrected.html')
        stats = create_3d_visualization_with_rna(
            cell_labels, nuclei_mask, rna_coords, cell_id, output_path,
            cell_smooth_sigma=args.cell_smooth,
            nucleus_smooth_sigma=args.nucleus_smooth,
            rna_sphere_size=args.rna_size
        )
        
        if stats:
            all_stats[cell_id] = stats
    
    # 保存统计信息
    stats_file = os.path.join(args.output, 'rna_statistics_corrected.json')
    with open(stats_file, 'w') as f:
        json.dump(all_stats, f, indent=2)
    print(f"\n✅ 统计信息已保存: {stats_file}")
    
    print("\n" + "=" * 60)
    print(f"✅ 完成！共创建 {len(all_stats)} 个 3D 可视化")
    print("\n修正后的 RNA 统计:")
    for cell_id, stats in all_stats.items():
        print(f"\n细胞 {cell_id}:")
        print(f"  原始 RNA: {stats['original_rna_count']}")
        print(f"  过滤后 RNA: {stats['filtered_rna_count']} (保留 {stats['filtered_rna_count']/stats['original_rna_count']*100:.1f}%)")
        print(f"  细胞核内: {stats['rna_in_nucleus']} ({stats['rna_in_nucleus']/max(stats['filtered_rna_count'], 1)*100:.1f}%)")
        print(f"  细胞质内: {stats['rna_in_cytoplasm']} ({stats['rna_in_cytoplasm']/max(stats['filtered_rna_count'], 1)*100:.1f}%)")
    print("=" * 60)


if __name__ == '__main__':
    main()
