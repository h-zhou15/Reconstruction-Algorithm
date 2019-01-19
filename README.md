# Reconstruction-Algorithm
## 图像重建学习中的算法code记录
## BP.m
  第一次学习解析重建算法，自己手动画了一个椭圆，模拟radon变换解了投影值方程<br>
  就是一个很简单的弦长公式。然后分别用了<br>
  直接反投影（效果最差）、<br>
  滤波反投影（窗函数为斜坡函数+矩形窗，效果更好，但是存在一点瑕疵，有空的时候改一下）<br>
  卷积反投影（斜坡函数在时域中的表现形式，直接卷积滤波，滤波函数在时域中较为复杂，不建议）<br>
  
## fbp.m

  核医学仪器与方法SPECT大作业<br>
  未对投影数据做修正，在忽略康普顿散射的情况下，发射成像可以近似等效于透射成像（SPECT、PET~=平行束CT）<br>
  原始文件为人体的全身骨投影数据，把三维投影数据分成128张二维投影切片做fbp<br>
  fbp选用HammingFUnction+斜坡函数做滤波器，时间很快，效果还行（对比ART）<br>
  
## Siddon.m
  迭代重建算法，采用线模型计算传输矩阵。其中快速计算相交线长度的code。<br>
  根据paper<br>
  Robert L. Siddon,1985,AAPM,
  Fast calculation of the exact radiological path for a three-dimensional CT array.<br>
  文章提出的针对三维CT数据快速计算相交线长度。我在这里简化为二维情况<br>
  根绝计算得到的相交线长度，生成了OSEM中需要的传输矩阵。<br>
  
## OSEM.m
  在siddon算法的基础上，以单线模型为基准的OSEM算法<br>
  应该重点理解一下EMML的迭代公式<br>
  https://blog.csdn.net/budong_2017/article/details/79779298
