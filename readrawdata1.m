%���ļ���ȡ׶��Բ���ƽ��̽����CT����
%����ԭʼ�ļ�Ŀ¼������
%ͶӰ���ݱ�����256*256*360����ά����p��

a=fopen('Circular CBCT_flat_panel_detector.prj');
p=fread(a,'float');
p=reshape(p,256,256,360);
fclose(a);