    def map_to_genomic_manually(self, db):
        # TODO: Not congruent w/ automatic for NM_015015.2:c.2441+1G>A
        name = self.tx
        c_pos = self.base
        offset = self.offset

        coding = []
        for i in db.children(f'rna-{name}', featuretype='CDS', order_by='start'):
            coding.append(i)
    
        strand = db[f'rna-{name}'].strand
        l = 0
        if strand == '-':
            coding = coding[::-1]  # reverse CDS order
            
            for i in coding:
                if l + len(i) >= c_pos:
                    break
                else:
                    l += len(i)
            # i.end holds genomic coordinate
            g_pos = i.end - (c_pos - l) + 1
        
        else:
            for i in coding:
                if l + len(i) >= c_pos:
                    break
                else:
                    l += len(i)
            g_pos = i.start + c_pos - 1

        if offset:
            # The offset is in relation to the coding sequence, eg "+3" means
            # three bases downstream of the previous coding sequence.
            if strand == '-':
                g_pos -= offset
            else:
                g_pos += offset

        # click.echo(log(f'Variant on chromosome {chromosome}, g. position {g_pos}'))
        return g_pos