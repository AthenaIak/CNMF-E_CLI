function[movieID, numFiles] = getMovieID(mouseID, session, day, trial, isRest)

stamps = {
    {32363, 1, 1, 1, false, '20170710_132114', 5};
};

selectedIdx = 1; %TODO :find it 

movieID = stamps{selectedIdx}{6};
numFiles = stamps{selectedIdx}{7};
