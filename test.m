A_1_cell = cell(1,3);                                                                                            % Source Strengths for the distributed point sources at the top of interface
for indice=1:3
    for i=1:NumSourcePt_IntrFc_y
        for j=1:NumSourcePt_IntrFc_x
            index=j+(i-1)*NumSourcePt_IntrFc_x+(indice-1)*NumSourceTot;
            A_1_cell{indice}(i,j)=A_1(index);
        end
    end
end