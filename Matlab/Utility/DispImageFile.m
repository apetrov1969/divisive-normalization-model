function [retval_handle,...
          retval_loadedImg] = DispImageFile(image_path)

%DispImageFile - Load and show an image
% 1.0.0 2014-02-18 TS -- Wrote it from scratch

retval_handle = figure('Name',image_path);
retval_loadedImg = imshow(image_path);

end %%% of file