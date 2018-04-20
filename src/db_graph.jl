using Bio.Seq
using Bio

function twin{k}(::Type{DNAKmer{k}}, km)
    DNAKmer{k}(Bio.Seq.reverse_complement(DNASequence("$km")))
end

function build{k}(::Type{DNAKmer{k}}, fastafiles)
  kmers = Set{DNAKmer{k}}()
  record = FASTASeqRecord{Bio.Seq.DNASequence}()
  for fastafile in fastafiles
    open(FASTAReader{Bio.Seq.DNASequence}, fastafile) do reader
      while !eof(reader)
        read!(reader, record)
        for (pos, kmer) in each(DNAKmer{k}, record.seq)
          push!(kmers, kmer)
        end
        for (pos, kmer) in each(DNAKmer{k}, Bio.Seq.reverse_complement(record.seq))
          push!(kmers, kmer)
        end
      end
    end
  end
  return kmers
end

function contig_to_string(c)
  return "$(c[1])" * join(["$(x[end])" for x in c[2:end]])
end

function get_contig{k}(::Type{DNAKmer{k}}, kmers, km)
  contig_fw = get_contig_forward(DNAKmer{k}, kmers, km)
  contig_bw = get_contig_forward(DNAKmer{k}, kmers, twin(DNAKmer{k}, km))

  if km in neighbors(contig_fw[end])
    return contig_fw
  else
    return [[twin(DNAKmer{k}, x) for x in reverse(contig_bw)[1:end-1]] ; contig_fw]
  end
end

function get_contig_forward{k}(::Type{DNAKmer{k}}, kmers, km)
    c_fw = DNAKmer{k}[km]
    while true
        # check if forward exists, and only 1:
        fw_neighbors = [kmer for kmer in neighbors(c_fw[end]) if kmer in kmers]
        if length(fw_neighbors) != 1
          break
        end
        candidate = fw_neighbors[1]
        if candidate == km || candidate == twin(DNAKmer{k}, km) || candidate == twin(DNAKmer{k}, c_fw[end])
          break
        end
        bw_neighbors = [kmer for kmer in neighbors(twin(DNAKmer{k}, candidate)) if kmer in kmers]
        if length(bw_neighbors) != 1
          break
        end
        push!(c_fw, candidate)
    end
    return c_fw
end

function all_contigs{k}(::Type{DNAKmer{k}}, kmers)
  done = Set{DNAKmer{k}}()
  contigs = []
  for kmer in kmers
    if kmer in done
      continue
    end
    contig = get_contig(DNAKmer{k}, kmers, kmer)
    for c_kmer in contig
      push!(done, c_kmer)
      push!(done, twin(DNAKmer{k}, c_kmer))
    end
    push!(contigs, contig)
  end
  return contigs
end

function db_graph_contig_kmers{k}(::Type{DNAKmer{k}}, fastafiles)
  kmers = build(DNAKmer{k}, fastafiles)
  # contig_kmers = Set{DNAKmer{k}}()
  contig_kmers = Dict{DNAKmer{k},Int16}()
  for contig in all_contigs(DNAKmer{k}, kmers)
      # push!(contig_kmers, canonical(contig[1]))
      contig_kmers[canonical(contig[1])] = length(contig)
  end
  # println(maximum(values(contig_kmers)))
  return contig_kmers
end


# # main:
#   k = 31
#   contigs = db_graph_contigs(DNAKmer{k}, ["/projects/pathogist/cgMLST/MTB/Rv0900.fa"])
#   for contig in contigs
#     println(contig_to_string(contig))
#   end
