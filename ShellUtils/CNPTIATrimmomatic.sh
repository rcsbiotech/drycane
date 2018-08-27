
# Pipeline de pré-processamento com Trimmomatic
# Os arquivos devem estar formatados da seguinte forma:
# SAMPLE01_R1.fastq			-> Arquivo forward
# SAMPLE02_R2.fastq			-> Sequência reverse
# Os nomes anteriores a "_R" serão mantidos na saída.

# Variáveis
# input 	= diretório com os arquivos no formato .fastq
# outdir 	= diretório de saída dos arquivos do Trimmomatic
# threads 	= número de processadores
# base_out 	= diretório de saída

input="$1"

if [ ! ${input} ]
then   
	echo "Missing input path to directory containing *.fastq files"
        exit
else   
        if [ ! -d ${input} ]
        then   
		echo "Wrong input path (${input}). Please, check the first argument."
                exit
        fi
fi

outdir="$2"

if [ ! ${outdir} ]
then   
	outdir="."
else
	if [ ! -d ${outdir} ]
	then
		echo "Output directory (${outdir}) doesn't exist. Please, choose another one."
		exit
	fi	
fi

threads="$3"

# criação dos diretórios de saída
# diretório base de saída:
base_out="${outdir}"

# diretório para os logs de saída:
log_dir="${base_out}/logs"

# diretório para os sumários:
summ_dir="${base_out}/summaries"

# diretório para depósito das sequências processadas:
proc_dir="${base_out}/output"

# criação dos diretórios de saídas:

mkdir -p ${log_dir}
mkdir -p ${summ_dir}
mkdir -p ${proc_dir}

# Listar todos os arquivos fastq (dgpinheiro)

for fastq1 in `ls ${input}/*_R1*.fastq`; do
	fastq2=`echo ${fastq1} | sed 's/R1/R2/'`;
	fqname1=`basename ${fastq1} .fastq`
	fqname2=`basename ${fastq2} .fastq`
	fqbase=`basename ${fastq1} .fastq | sed 's/_R1//'`
	
	java -jar /home/externo/danielgp/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
	-threads ${threads} \
	${fastq1} \
	${fastq2} \
	-trimlog ${log_dir}/${fqbase}.log_file.txt \
	-summary ${summ_dir}/${fqbase}.summary_file.txt \
	-baseout ${proc_dir}/${fqbase} \
	ILLUMINACLIP:/home/externo/danielgp/adapters.fa:2:40:15 \
	LEADING:5 \
	TRAILING:5 \
	SLIDINGWINDOW:4:30 \
	MINLEN:35
	
	done
	
