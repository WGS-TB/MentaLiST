using Bio.Seq
using Bio

function twin{k}(::Type{DNAKmer{k}}, km)
    DNAKmer{k}(Bio.Seq.reverse_complement(DNASequence("$km")))
end

function build{k}(::Type{DNAKmer{k}}, fastafiles)
  kmers = Set{DNAKmer{k}}()
  record = FASTASeqRecord{BioSequence{DNAAlphabet{2}}}()
  for fastafile in fastafiles
    open(FASTAReader{BioSequence{DNAAlphabet{2}}}, fastafile) do reader
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
  contig_kmers = Dict{DNAKmer{k},Int16}()
  for contig in all_contigs(DNAKmer{k}, kmers)
      contig_kmers[canonical(contig[1])] = length(contig)
  end
  return contig_kmers
end


# Calling:
function sequence_coverage{k}(::Type{DNAKmer{k}}, sequence, kmer_count, kmer_thr=6)
  # find the minimum coverage depth of all kmers of the given variant, and
  # how many kmers are covered/uncovered
  smallest_coverage = 100000
  covered = 0
  uncovered = 0
  for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
    can_kmer = canonical(kmer)
    cov = get(kmer_count, can_kmer, 0)
    if cov >= kmer_thr
      covered += 1
    else
      uncovered += 1
    end
    if cov < smallest_coverage
      smallest_coverage = cov
    end
  end
  return smallest_coverage, covered, uncovered
end

function find_allele_variants{k}(::Type{DNAKmer{k}}, allele_seqs, template_alleles, kmer_count, kmer_thr=8, max_mutations=10)
  novel_allele_list = []
  found = Set{String}()
  function scan_variant(n_mut, allele, sequence, events, start_pos=1)
    if n_mut > max_mutations
      return
    end
    last_kmer_present = false
    for (pos, kmer) in each(DNAKmer{k}, DNASequence(sequence))
      if pos < start_pos
        continue
      end
      can_kmer = canonical(kmer)
      if get(kmer_count, can_kmer, 0) >= kmer_thr
        if pos > start_pos +1 && (!last_kmer_present)
          # found a kmer and the last was not there; check if I can find kmer variations for the previous position
          km = "$kmer"
          for base in ["A","C","T","G"]
            if get(kmer_count, canonical(DNAKmer{k}("$base$(km[1:end-1])")), 0) >= kmer_thr
              # found a kmer variant;
              # Test one more kmer before adding to the pool:
              # substitution at position pos-1
              substitution_seq = sequence[1:pos-2] * base * sequence[pos:end]
              push!(candidate_list, (n_mut+1, allele, substitution_seq, [events; [("S",pos-1,"$(sequence[pos-1])->$base")]], start_pos))
              # insertion:
              insertion_seq = sequence[1:pos-1] * base * sequence[pos:end]
              push!(candidate_list, (n_mut+1, allele, insertion_seq, [events; [("I",pos-1,base)]], start_pos))
              # deletion:
              for i in 0:2
                if pos-2-i < 1
                  break
                end
                if sequence[pos-2-i:pos-2-i] == base
                  deletion_seq = sequence[1:pos-2-i] * sequence[pos:end]
                  push!(candidate_list, (n_mut+1+i, allele, deletion_seq, [events; [("D",pos-1,"$(i+1)")]], start_pos))
                  break
                end
              end
              # this variant has mismatches, so just return;
              return
            end
          end
        end
        last_kmer_present = true
      else
        if last_kmer_present
          # last_kmer_present was present, but this not; check kmer variations for nest position
          km = "$kmer"
          for base in ["A","C","T","G"]
            if get(kmer_count, canonical(DNAKmer{k}("$(km[1:end-1])$base")), 0) >= kmer_thr
              # println("Found a $base at pos $(pos+k-1)")
              # substitution at position pos+k-1;
              substitution_seq = sequence[1:pos+k-2] * base * sequence[pos+k:end]
              push!(candidate_list, (n_mut+1, allele, substitution_seq, [events; [("S",pos+k-1,"$(sequence[pos+k-1])->$base")]], pos-1))
              # insertion:
              insertion_seq = sequence[1:pos+k-2] * base * sequence[pos+k-1:end]
              push!(candidate_list, (n_mut+1, allele, insertion_seq, [events; [("I",pos+k-1,base)]], pos-1))
              # deletion:
              for i in 0:2
                # println("trying $(pos+k+i), found $(sequence[pos+k+i]) (want $base) $(sequence[pos+k+i:pos+k+i] == base).")
                if pos+k+i > length(sequence)
                  break
                end
                if sequence[pos+k+i:pos+k+i] == base
                  # println("Found a possible deletion, pos $(pos+k-1) to $(pos+k+i-1)")
                  deletion_seq = sequence[1:pos+k-2] * sequence[pos+k+i:end]
                  push!(candidate_list, (n_mut+1+i, allele, deletion_seq, [events; [("D",pos+k-1,"$(i+1)")]], pos-1))
                  break
                end
              end
              # the original sequence has mismatches, so just return;
              return
            end
          end
        end
        last_kmer_present = false
      end
    end
    # passed without suggesting mutations, test this variant:
    var_ab = sequence_coverage(DNAKmer{k}, sequence, kmer_count)[1]
    if var_ab >= kmer_thr
      if !in(sequence, found)
        push!(novel_allele_list, (n_mut, allele, sequence, events, var_ab))
        push!(found, sequence)
        max_mutations = n_mut
      end
    end
  end
  # candidate_list = [(0, template_sequence, [], 0) for template_sequence in template_sequences]
  candidate_list = [(0, allele, allele_seqs[allele], [], 0) for allele in template_alleles]
  idx = 1
  while idx <= length(candidate_list)
    n_mut, allele, candidate_sequence, mutations, start_pos = candidate_list[idx]
    scan_variant(n_mut, allele, candidate_sequence, mutations, start_pos)
    idx += 1
  end

  # println("Variants found: $(length(novel_allele_list))\n") #$novel_allele_list")
  # println("Variants found: $(length(novel_allele_list))\n$novel_allele_list")
  if length(novel_allele_list) > 0
    # sort and return the one with minimum mutations: TODO: tie?
    return sort(novel_allele_list, by=x->x[1])[1]
  else
    return nothing
  end

end



function describe_mutation(mut)
  mut, pos, desc = mut
  if mut == "D"
    return "Del of len $desc at pos $pos"
  elseif mut == "I"
    return "Ins of base $desc at pos $pos"
  elseif mut == "S"
    return "Subst $desc at pos $pos"
  else
    return ""
  end
end
