import json

r = json.load(open("/home/manec/vcc-avaliador/peers-k4.json"))
for v in r["variantes"]:
    n = v["variante"]
    if "aborta_antes_do_run_meta" in v:
        print(f"{n:24s} ABORTA: {v['aborta_antes_do_run_meta'][:170]}")
        continue
    print(f"{n:24s} oito={v['veredito_check_submission']:9s} "
          f"extra={v['divergem_dos_portoes_extra']} "
          f"ancora={v['veredito_expect_from_run_meta'][:110]}")
print()
print("valores de de_real_fornecida:", r["variantes"][-2]["valores"])
