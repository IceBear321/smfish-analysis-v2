#!/usr/bin/env python3
"""
从 CZI 文件中提取所有通道并保存为单独的 TIFF 文件
"""

import argparse
import os
import numpy as np
import czifile
from skimage import io
import matplotlib.pyplot as plt


def extract_channels(czi_path, output_dir):
    """从 CZI 文件提取所有通道"""
    
    print(f"读取 CZI 文件: {czi_path}")
    with czifile.CziFile(czi_path) as czi:
        # 读取图像数据
        image_data = czi.asarray()
        print(f"  原始形状: {image_data.shape}")
        
        # 通常 CZI 格式为 (T, Z, C, Y, X) 或类似
        # 我们需要找到 C (通道) 维度
        if image_data.ndim == 5:
            # (T, Z, C, Y, X)
            image_data = image_data[0]  # 取第一个时间点
            print(f"  提取后形状: {image_data.shape}")
        
        if image_data.ndim == 4:
            # (Z, C, Y, X)
            num_channels = image_data.shape[1]
            print(f"  检测到 {num_channels} 个通道")
            
            os.makedirs(output_dir, exist_ok=True)
            
            # 创建预览图
            fig, axes = plt.subplots(1, num_channels, figsize=(5*num_channels, 5))
            if num_channels == 1:
                axes = [axes]
            
            # 提取每个通道
            for c in range(num_channels):
                channel_data = image_data[:, c, :, :]
                
                # 保存为 TIFF
                output_path = os.path.join(output_dir, f'channel_{c}.tif')
                io.imsave(output_path, channel_data.astype(np.uint16), check_contrast=False)
                print(f"  ✓ 通道 {c} 已保存: {output_path}")
                print(f"    形状: {channel_data.shape}")
                print(f"    数据类型: {channel_data.dtype}")
                print(f"    范围: {channel_data.min()} - {channel_data.max()}")
                print(f"    均值: {channel_data.mean():.2f}")
                
                # 添加到预览图 (MIP)
                mip = channel_data.max(axis=0)
                axes[c].imshow(mip, cmap='gray')
                axes[c].set_title(f'Channel {c}\n(mean={channel_data.mean():.2f})')
                axes[c].axis('off')
            
            # 保存预览图
            preview_path = os.path.join(output_dir, 'channels_preview.png')
            plt.tight_layout()
            plt.savefig(preview_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"\n✓ 预览图已保存: {preview_path}")
            
        else:
            print(f"  ⚠️  不支持的图像维度: {image_data.ndim}")
            return
        
        # 提取元数据
        try:
            metadata = czi.metadata()
            
            # 尝试提取像素尺寸
            import xml.etree.ElementTree as ET
            root = ET.fromstring(metadata)
            
            # 查找 Scaling 元素
            scaling = root.find('.//Scaling')
            if scaling is not None:
                items = scaling.findall('.//Distance')
                pixel_sizes = {}
                for item in items:
                    dim_id = item.get('Id')
                    value = item.find('Value')
                    if value is not None:
                        pixel_sizes[dim_id] = float(value.text) * 1e6  # 转换为微米
                
                print(f"\n像素尺寸:")
                for dim, size in pixel_sizes.items():
                    print(f"  {dim}: {size:.6f} μm/pixel")
                
                # 保存像素尺寸
                pixel_size_file = os.path.join(output_dir, 'pixel_sizes.txt')
                with open(pixel_size_file, 'w') as f:
                    for dim, size in pixel_sizes.items():
                        f.write(f"{dim}: {size:.6f} μm/pixel\n")
                print(f"✓ 像素尺寸已保存: {pixel_size_file}")
        
        except Exception as e:
            print(f"⚠️  无法提取元数据: {e}")
    
    print(f"\n✅ 完成！所有通道已提取到: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='从 CZI 文件中提取所有通道'
    )
    parser.add_argument('input', help='输入 CZI 文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    
    args = parser.parse_args()
    
    extract_channels(args.input, args.output)


if __name__ == '__main__':
    main()
