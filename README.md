# Spectra analysis with eigendecomposition

Set of spectra can be organised column-wise into an absorbance matrix, **W**. From now on, each column of W shall be treated as </br>a vector **w**<sub>i</sub> in the m-dimensional vector space. According to the Lambert-Beer law, the absorption is linearly dependent on concentration. This means that the absorbance matrix **W** can be represented as a product of the two matrices: the concentration matrix **C** and the matrix of basic spectra **B**. In order to decompose matrix **W** into matrices that contain the concentration and spectra of the pure components respectively, the absorbance matrix must be verified and properly prepared. Verifying the data includes analysis of: mean and standard deviation values for each column, as well as correlation coefficients analysis of each  pair of columns. The curve showing the mean values and standard deviations of each vector should be reasonably smooth and there should be no sharp spikes in the correlation coefficients plot. All spectra (vectors) that produce a spike on one of the curves mentioned above should be removed from the data set.

## Eigenvector decomposition
The covariance matrix W is decomposed into singular values (according to the SVD algorithm), which allowed to determine eigenvalues and eigenvectors (U). The contributions of the principal components were calculated using the matrix 
equation:

 </br>$U = [V_{T}V]^{-1}V^{T}W$ </br>
 
 Conversely, the coordinates of the individuals in the space of principal components are 
calculated by:

 </br>$V = WU[UU^{T}]^{-1}$

V – the coordinates of the individuals in the space of principal components,  
W - centered spectrum,  
U – contributions of the principal components.
</br>Matrices U and V were used to calculate the components of the spectra. If the matrix 
of the components of the spectra is denoted by C, it can be calculated as: 

</br>$C = VU^{T}$ 

</br>Based on C, reconstructed spectra (D) were determined and then residual spectra (R) 
were calculated by: 

</br>$R = W - D$ 

</br>The spectra are placed in a vector space, which is spanned by the principal 
components. These vectors are placed into a hyperplane with the number
of dimensions equal to the number of pure spectral forms. This means that the first 
step is to determine the number of spectral forms. It can be done by the analysis of the 
eigenvalues and residual spectra. It is assumed that the rank of residual spectra
at which the systematic changes of these spectra disappear defines the number
of significant principal components (spectral forms). Subsequently, the points 
corresponding to the pure spectral forms need to be found over the principal 
component space.
## Possible options
This code allows a preliminary analysis of the quality of the recorded spectra. It allows the evaluation of UV-Vis, CD and IR spectra. You can also verify the mean values, standard deviations and correlation coefficients between spectra. It also allows you to view centered and standardized spectra, as well as components, reconstructed spectra and residual spectra resulting from spectra eigendecomposition.
## How to run the program
Simply run the program with the python3 command.
