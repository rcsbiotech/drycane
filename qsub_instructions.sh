
# How to: QSUB

echo "command" | qsub \
-cwd \
-N "name" \
-m base -M rcs.biotec@gmail.com \
mf=30G

# example usage

echo "bash /home/externo/danielgp/bin/rsem-eval-rcsilva.sh" | qsub \
-N alignOK -cwd -l mf=30G


