from pathlib import Path
p=Path(__file__).parent/'stage/lapw_valence_state.py';s=p.read_text();needle='''        raise ValueError("assembled raw projectors differ from native snapshot")''';assert needle in s;s=s.replace(needle,'''        import os
        audit = os.environ.get("EDMFT_RESPONSE_RAW_AUDIT")
        if audit:
            np.savez(audit + "_" + case + spin_channel + ".npz",
                assembled=np.array([r[:nret] for r in raw]), expected=model.raw,
                energies=model.eps)
        difference=np.linalg.norm(np.array([r[:nret] for r in raw])-model.raw)/np.linalg.norm(model.raw)
        raise ValueError(f"assembled raw projectors differ from native snapshot: relative error {difference:g}")''');p.write_text(s)
