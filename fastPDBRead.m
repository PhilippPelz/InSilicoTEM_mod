function [PDBdata] = fastPDBRead(readFile)
%%  -- Manual --
% This program is the most speedy way to read a PDB file that I could come
% up with. It's function is simple: give it a PDB file and out comes a
% matlab-friendly set of data. In cumbersomely large PDB's (such as those that 
% include solvent), this can shave off a good amount of time relative to
% many programs. Unfortunately there is no easy way to hasten the slowest
% step, which is turning strings into doubles.
%
% It outputs 10 pieces of information about the PDB:
% 
% atomType (the class or type of atom, such as ATOM, HETATM, SOL, etc)
% atomNum  (index number of the atom)
% atomName (elemental identification of the atom)
% resName  (name of the amino acid/residue)
% resNum   (index number of the amino acid)
% 
% X        (X position of atom)
% Y        (Y position of atom)
% Z        (Z position of atom)
% b_factor (this is extra information about an atom. In PDBQT's it holds 
%           the partial charge, for CHARMM this is the chain name, and so on)
% comment  (anything left on each line. In charmm this is the chain name,
%           in many other programs this is the element and charge)
%
% example usage: plot the atoms of 3IJU.pdb
% 
% atoms = fastPDBRead('3IJU.pdb')
% plot3(atoms.X, atoms.Y, atoms.Z, '.');
% 
%% --- HOW TO MAKE THIS CODE FASTER! --- >> COMMENT OUT WHAT YOU DON'T USE!!
%
% This program reads everything about the PDB by default. If you want a
% faster code for whatever reason, you can comment out the lines you don't
% need. Each numeric data removed (such at resNum, or b_factor) speeds it 
% up by 7-8%. Each string data removed (such as resName or atomName) speeds
% it up by 1-2%.

%% -- OUTPUT --

tic;

% initialize file
FileID = fopen(readFile);
rawText = fread(FileID,inf,'*char');

% parse lines by end-of-lines
splitLines = strread(rawText, '%s', 'delimiter', '\n');

% initialize variables
numLines = length(splitLines);

PDBdata.atomType = cell(1,numLines);
PDBdata.atomNum  = cell(1,numLines);
PDBdata.atomName = cell(1,numLines);
PDBdata.resName  = cell(1,numLines);
PDBdata.resNum   = cell(1,numLines);
PDBdata.X        = cell(1,numLines);
PDBdata.Y        = cell(1,numLines);
PDBdata.Z        = cell(1,numLines);
PDBdata.b_factor = cell(1,numLines);
PDBdata.comment  = cell(1,numLines);
PDBdata.element  = cell(1,numLines);
PDBdata.occupancy  = cell(1,numLines);

% read each line
m = 1;
for n = 1:numLines
    
    thisLine = cell2mat(splitLines(n));
%     fprintf('%s\n',thisLine(1:4))
    if length(thisLine) > 60 && sum(isstrprop(thisLine(23:66), 'alpha')) == 0 && strcmp(thisLine(1:4),'ATOM')
        
        PDBdata.atomType(m) = {thisLine(1:6)};
        PDBdata.atomNum(m)  = {thisLine(7:11)};
        PDBdata.atomName(m) = {thisLine(13:16)};
        PDBdata.resName(m)  = {thisLine(18:20)};
        PDBdata.resNum(m)   = {thisLine(23:26)};
        PDBdata.X(m)        = {thisLine(31:38)};
        PDBdata.Y(m)        = {thisLine(39:46)};
        PDBdata.Z(m)        = {thisLine(47:54)};
        PDBdata.b_factor(m) = {thisLine(61:66)};
        PDBdata.comment(m)  = {thisLine(67:end)};
        PDBdata.element(m)  = {thisLine(78:79)};
        PDBdata.occupancy(m)  = {thisLine(57:60)};
        m = m + 1;
    end
    
end

% trim exess
PDBdata.atomType(m:end) = [];
PDBdata.atomNum(m:end)  = [];
PDBdata.atomName(m:end) = [];
PDBdata.resName(m:end)  = [];
PDBdata.resNum(m:end)   = [];
PDBdata.X(m:end)        = [];
PDBdata.Y(m:end)        = [];
PDBdata.Z(m:end)        = [];
PDBdata.b_factor(m:end) = [];
PDBdata.comment(m:end)  = [];

% reformat data for convenience
PDBdata.atomType = strtrim(PDBdata.atomType);
PDBdata.atomNum  = str2double(PDBdata.atomNum);
PDBdata.atomName = strtrim(PDBdata.atomName);
%PDBdata.resName  = strtrim(PDBdata.resName); % never needed in valid PDB
PDBdata.resNum   = str2double(PDBdata.resNum);
PDBdata.X        = str2double(PDBdata.X);
PDBdata.Y        = str2double(PDBdata.Y);
PDBdata.Z        = str2double(PDBdata.Z);
PDBdata.b_factor = str2double(PDBdata.b_factor);
PDBdata.comment  = strtrim(PDBdata.comment);
PDBdata.element  = strtrim(PDBdata.element);
PDBdata.occupancy = str2double(PDBdata.occupancy);
% close file
fclose(FileID);

toc;

end