function fixpermil(oldfile,newfile)

% FIXPERMIL  A function to fix notation for permil in postscript Matlab figures
%  
%   fixpermil(oldfile,newfile)
%
% Matlab improperly codes the symbol for \permil in some postscript.  This function
% edits the raw postscript automatically so that the the symbol, called as char(8240), 
% will display correctly in the Postscript file.
% 
% This function should be called after a 'print' command which has used the -depsc2 switch
%
% Details on this MATLAB bug here:
%
% http://www.mathworks.com/matlabcentral/answers/1231-does-matlab-2010b-intentionally-omit-char-8240-and-others-or-is-this-a-bug-with-a-fix
%
% Note that on my system, the incorrect \permil coding is \32, and needs to
% be replaced with \275.  It is possible this is version specific. 
%
% Example:
% figure(1); clf; plot(randn(100,1),randn(100,1)); ylabel(['5 ',char(8240)]);
% print -depsc2 figure.eps
% fixpermil('figure.eps','newfigure.eps')
%
% Kevin Anchukaitis, kja@whoi.edu
%
% Revision history: 
% First version 07.12.08, revised 12.09.09
% extensively rewritten 09.23.13 to use strrep, based on improvements by Patrick Sturm at ETH
% http://homepage.agrl.ethz.ch/~pasturm/eddycalc/otherMATLABfunctions/html/fixPermil.html

% open the file and read it in as a long string
fid = fopen(oldfile,'r');
str = fread(fid);
str = char(str');

% if only a single argument was passed, use the old filename as a the new filename
if nargin ==1; newfile=oldfile; end

% First, deal with the ISOLatinEncoding, new based on Patrick Sturm's modifications
[id1] = findstr(str,'/ISOLatin1Encoding');
for i=1:length(id1)
  str = strrep(str,'/ISOLatin1Encoding','/StandardEncoding');
end
           
% Then, replace \32 with \275 to give the proper encoding for permil
% It is possible this encoding changes over operating systems of MATLAB versions
[id2]   = findstr(str, '\32');
for i = 1:length(id2)
str = strrep(str,'\32','\275');
end

% write out the new file
fid = fopen(newfile, 'w');
fprintf(fid, '%s', sprintf('%s',str));
fclose(fid);
