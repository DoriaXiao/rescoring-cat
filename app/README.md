# CAT Rescoring Explorer (companion app)

Interactive Gradio app accompanying the paper. Lets readers explore the
simulation under user-specified parameters across four tabs:

1. **Cross-Population Misalignment** — scoring vs selection error across ρ
2. **Random Noise & Shifts** — random calibration noise and structured difficulty shifts
3. **Stepwise: Short Tests Benefit Most** — the rescoring advantage by test length
4. **Single-Point Explorer** — any combination of ρ and test length

**Live app:** https://huggingface.co/spaces/doriax/cat-rescoring-explorer

## Run locally

```bash
pip install gradio numpy matplotlib
python app.py
```

## Deploy on Hugging Face Spaces

Create a new Space with the **Gradio** SDK and add `app.py` plus a
`requirements.txt` containing:

```
numpy
matplotlib
gradio>=4.0
```

The app re-implements the IRT/CAT simulation in NumPy (no R dependency), so it
runs standalone. It uses much smaller simulations than the paper for speed
(the paper uses 100 examinees × 50 replications per condition), so estimates at
the default settings are noisy.
