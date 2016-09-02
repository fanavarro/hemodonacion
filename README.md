# hemodonacion
El presente repositorio contiene los scripts y los datos obtenidos para el trabajo final del máster de bioinformática de la Universidad de Murcia.

El repositorio consta de las siguientes carpetas:
	src -> contiene los scripts en Perl y R desarrollados.
		ensemblMining.pl es el script de Perl encargado de obtener un fichero csv (cuyo nombre se especifica por parametro) con todas las mutaciones que afectan al codon iniciador de la traduccion. 
		statistics.R es el script de R usado para analizar el fichero csv obtenido mediante el script ensemblMining.pl.
		exomes.R es el script de R usado para analizar los exomas de 7 pacientes individuales. Deben ejecutarse las instrucciones de statistics.R para que quede en memoria el objeto "csv", usado en este script.


	data/exomes -> contiene ficheros de tipo xlsx con la información de los exomas de 7 pacientes individuales.
	data/tsv -> contiene el fichero "final_out_no_filter.csv", con el fichero obtenido mediante el script ensemblMining.pl y el fichero "final_out_filter.csv", resultado de filtrar el primero para eliminar transcritos que no producen proteina.

Para ejecutar ensemblMining.pl es necesario:
	Perl 5.10
	API de ensembl para perl.
	Librerias de Perl:
		LWP::Simple::REST
		LWP::Simple
		JSON
		XML::LibXML
		HTML::TableExtract
		Class::Singleton
		List::Util
