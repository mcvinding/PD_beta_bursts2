function [thickavg, thickstd] = read_fs_thickness(fname)
% Read the ROUJ cortical thickness from Freesurfer output stat file.
% USE: [thickavg, thichstd] = read_fs_thickness(fname)

fid = fopen(fname,'r');
linenum = 53;
C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
c = (C{1,1});
c = c{1};
X = splitstring(c);

% fid = fopen(fname,'r');
% linenum = 52;
% H = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
% h = char(H{1,1})
% Y = splitstring(h)
% Y = Y(3:end,:)

thickavg = str2double(X(5,:));
thickstd = str2double(X(6,:));

fprintf('Average thickness %.4f (sd %.4f)\n', thickavg, thickstd)
end