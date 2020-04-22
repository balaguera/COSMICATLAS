// DOCUMENTATION FOR DOXYGEN

/**
 * @mainpage
 *
 * @section intro_sec Introduction
 COSMICATLAS is a cosmological C++ package with two main branches. One one side, Bias Assignment Method (BAM), aimed at generatting mock catalogs.<br>
On the other hand, the CosmicLib code, aimed at producing differnet theoretical curves of cosmological observables (e.g abundance, clustering). <br>




 * @subsection _author Authors
 <a href="https://abalant.wixsite.com/abalan"> Andrés Balaguera-Antolínez</a> & <a href = "https://franciscokitaura.wixsite.com/home">Francisco-Shu Kitaura</a>.<br>
 Instituto de Astrofísica de Canarias


 * @subsection _pres Presentations
Check <a href = "../Headers/TalkBam.pdf">here</a> for the most recent talk on BAM (Innsbrück 2020)
 
 * @subsection _colla Collaboratos
 Chia-Hsun Chuang (KIPAC)<br>
 Gustavo Yepes (UAM)<br>
 Marcos Pellejero (DIPC)<br>
 Cheng Zhao (EPFL)<br>
 Ariel Sánchez (MPE)<br>
 Martha Lippich (MPE)<br>
 K. Nagamine (Osaka U.)<br>
 Metin Ata (KIP-Tokyo)<br>
 Claudio Dalla Vecchia (IAC)<br>
 Raúl Angulo (DIPC)<br>
 

 
 * @subsection pap Publications:
The main references for BAM are <br>
<a href=" http://adsabs.harvard.edu/abs/2019MNRAS.483L..58B">Paper I: The bias Assigment Method</a><br>
<a href="http://adsabs.harvard.edu/abs/2019arXiv190606109B">Paper II: One simulation to Have them all</a><br>
 
<table>
<tr>
<td>
*@section ex Executing the Code 

* @subsection rr Running the code BAM
This file is compiled with several options
*@code
make clean;
make bam
./cosmicatlass.exe -option parameter_bam.ini
*@endcode

* @subsection rc Running the code Cosmiclib
* The user can also compile the cosmic_lib code
*@code
make clean;
make cosmiclib
./cosmiclib.exe parameter_cosmiclib.ini
*@endcode
</td>

<td>
<div class="column"></div>
<img src="../Headers/flow.jpg" alt="BAM Flowchart" width="600" height="700">
</td>

</tr>
</table>

*/
