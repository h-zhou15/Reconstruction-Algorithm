%本文件读取锥束圆轨道平板探测器CT数据
%放在原始文件目录下运行
%投影数据保存在256*256*360的三维矩阵p中

a=fopen('Circular CBCT_flat_panel_detector.prj');
p=fread(a,'float');
p=reshape(p,256,256,360);
fclose(a);