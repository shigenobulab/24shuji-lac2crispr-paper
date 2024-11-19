require 'optparse'

opt = OptionParser.new
bam, genome, target_pam_position = nil
opt.on('-b', '--bam BAM', "bam file"){|v|
  bam = v    # "BAM_Amp180802/S9.on.ApisBuc1.genome.bowtie2.sorted.bam")
}
opt.on('-g', '--genome FASTA', "genome in fasta format"){|v|
  genome = v #"../../../Data/Acyrthosiphon_pisum/NIBB_ApisBuc1/blastdb/ApisBuc1.genome.fa"
}
opt.on('-p', '--pam LOCATION', "PAM position in the following format: 'GL350588:97404-97406(-)'"){|v|
  target_pam_position = v #"GL350588:97404-97406(-)"
}

opt.parse!(ARGV)

#=== conf ===
#bam = ARGV[0]    # "BAM_Amp180802/S9.on.ApisBuc1.genome.bowtie2.sorted.bam")
#genome = ARGV[1] #"../../../Data/Acyrthosiphon_pisum/NIBB_ApisBuc1/blastdb/ApisBuc1.genome.fa"
#target_pam_position = ARGV[2] #"GL350588:97404-97406(-)"

#===

## Analyze target position
class CrisprPam
  def initialize(pam_position) # (ex) "GL350588:97404-97406(-)"
    m = /(.+):(\d+)-(\d+)\((.?)\)/.match(pam_position)
    target_scaffold, target_from, target_to, target_strand = m[1], m[2], m[3], m[4]
    target_cleavage = nil
    if target_strand == "+"
      target_cleavage = target_from.to_i - 3
    elsif target_strand == "-"
      target_cleavage = target_to.to_i + 3
    else
      raise
    end
    @cleavage_site = [target_scaffold, target_cleavage, target_strand]
  end
  
  attr_reader :cleavage_site

end


pileup_file = File.dirname(bam) + "_" + File.basename(bam) + ".plup"

pam = CrisprPam.new(target_pam_position)

region = "#{pam.cleavage_site[0]}:#{pam.cleavage_site[1] - 30}-#{pam.cleavage_site[1] + 30}"

cmd = "samtools mpileup -f #{genome} -x -Q 13 -A --no-BAQ -d 1000000 -r #{region} #{bam}  > #{pileup_file}"
STDERR.puts cmd
system cmd


#====

pileup_inspect_data =[]
pileup_inspect_file = pileup_file + ".inspect"
o = File.open(pileup_inspect_file, "w")

cols = %w{chr pos nreads match mismatch ins del %ins %del %indel detail}

o.puts "#" + cols.join("\t")

File.open(pileup_file).each do |l|
  a = l.chomp.split(/\t/)
  pup = a[4]
  str = pup.dup
  bases = []
  while str.size > 0
    if /[ACGTNacgtn\.\,]/.match(str[0])
      bases << str[0]
#      p bases
      str = str[1..-1]
    elsif str[0] == '*'
      bases << str[0]
      str = str[1..-1]
    elsif /\^.[ATGCNatgcn\.\,]/.match(str[0, 3])
#      p str[0,3]
      bases << str[0, 3][-1]
 #     p bases.last
      str = str[3..-1]
    elsif str[0] == '$'
      # skip
      str = str[1..-1]
    elsif m = /^\-([0-9]+)/.match(str)
      type = "del"
      len = m[1].to_i
      bases << str[0, len + 2]
      str = str[(len + 2)..-1]
    elsif m = /^\+([0-9]+)/.match(str)
      type = "ins"
      len = m[1].to_i
      bases << str[0, len + 2]
      str = str[(len + 2)..-1]

    else
      p str[0,10]
      raise
    end
  end
#  p bases.size
#  p bases
  counter = Hash.new(0)
  bases.each do |b|
    c = b[0].upcase
    counter[c] += 1
  end

  chr = a[0]
  pos = a[1]
  ref = a[2]
  nreads = a[3].to_i
  match = counter["."] + counter[","]
  mismatch = counter["A"] + counter["T"] + counter["G"] + counter["C"]
  ins = counter["+"]
  del = counter["*"] + counter["-"]

  indel = ins + del
#  perc_mismatch = mismatch.to_f / (nreads+1) * 100.0
  perc_ins = ins.to_f / (nreads) * 100
  perc_del = del.to_f / (nreads) * 100
  perc_indel = indel.to_f / (nreads) * 100.0

  detail =  [",", ".", "A", "T", "G", "C", "-", "*", "+"].map{|x| "#{x}:#{counter[x]}"}.join("|")

  data = [chr, pos, nreads , 
        match, mismatch, ins, del, 
        sprintf("%.2f", perc_ins),  sprintf("%.2f", perc_del),
        sprintf("%.2f", perc_indel),
        detail]
  o.puts data.join("\t")
  pileup_inspect_data << data
end
o.close

#====

indel_peak_file = pileup_file + ".indel_peak.txt"
o = File.open(indel_peak_file, "w")

## Analyze target
m = /(.+):(\d+)-(\d+)\((.?)\)/.match(target_pam_position)
target_scaffold, target_from, target_to, target_strand = m[1], m[2], m[3], m[4]
target_cleavage = nil
if target_strand == "+"
  target_cleavage = target_from.to_i - 3
elsif target_strand == "-"
  target_cleavage = target_to.to_i + 3
else
  raise
end

data_target_idx = pileup_inspect_data.find_index{|d| d[1] == target_cleavage.to_s}
tmpary = []
((data_target_idx - 10)..(data_target_idx + 10)).each  do |i|
  dat = pileup_inspect_data[i]
  nreads = dat[2].to_i
  ins = dat[5].to_i
  del = dat[6].to_i
  indel_perc = dat[9].to_f
  tmpary << [i, indel_perc, nreads, ins, del]
end

tmpary_sorted = tmpary.sort{|a, b| a[1] <=> b[1] }
dat_max_pos = pileup_inspect_data[tmpary_sorted[-1][0]]

p dat_max_pos

o.puts "#" + %w{bam_file %indel position distance_from_cleavage_site num_total_reads num_indel_reads}.join("\t")
out = [pileup_file, #File.basename(bam),
       dat_max_pos[9],
       "#{dat_max_pos[0]}:#{dat_max_pos[1]}",
       dat_max_pos[1].to_i - target_cleavage,
       dat_max_pos[2],
       dat_max_pos[5] + dat_max_pos[6]
]
o.puts out.join("\t")

puts "#" + %w{bam_file %indel position distance_from_cleavage_site num_total_reads num_indel_reads}.join("\t")
puts out.join("\t")

o.close

#=== plot indel percentage

#p pileup_inspect_data
rscript_file = pileup_file + ".indel_peak.plot.R"
rscript_imgout_png = pileup_file + ".indel_peak.plot.png"
rscript_imgout_pdf = pileup_file + ".indel_peak.plot.pdf"

rscript = %Q{
dat <- read.delim("#{pileup_inspect_file}")
png("#{rscript_imgout_png}")
plot(NA, type="n", ylim=c(0,100), xlim=c(0,62), ylab="%indels", main="#{pileup_inspect_file}")
lines(dat[10], type="o", ylim=c(0,100), col="royalblue", pch=16)
abline(v=31, lty=2, col="grey")
mtext(text=paste("peak: ",max(dat[10]), "%  "), adj=1, line=-2)
dev.off()

pdf("#{rscript_imgout_pdf}")
plot(NA, type="n", ylim=c(0,100), xlim=c(0,62), ylab="%indels", main="#{pileup_inspect_file}")
lines(dat[10], type="o", ylim=c(0,100), col="royalblue", pch=16)
abline(v=31, lty=2, col="grey")
mtext(text=paste("peak: ",max(dat[10]), "%  "), adj=1, line=-2)
dev.off()
}

File.open(rscript_file, "w"){|o| o.puts rscript}
cmd = "Rscript --slave --vanilla #{rscript_file} "
system cmd
