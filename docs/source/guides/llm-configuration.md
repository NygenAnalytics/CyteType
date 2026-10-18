---
description: Run CyteType's annotation agents on your own model provider, including AWS Bedrock, OpenAI-compatible endpoints and local Ollama models, and keep your provider keys out of files.
---

(llm-configuration)=
# Use your own LLM

By default, CyteType runs every job on models it hosts. With `llm_configs`, the annotation runs on models from your own provider account instead, for example to work under your institution's agreement with a provider, to use a particular model, or to use a model you host yourself.

What changes, and what does not:

- **Your models run all of the job's agents** (see [Several models](#llm-configuration-agents)). CyteType's hosted models are not mixed in for them.
- **CyteType's servers call your provider** during the job, with the key you pass. The calls count against your provider account.
- **Your [usage limits](../platform/limits.md) still apply,** to jobs with your own models as to any other job.
- **Chat and target analysis always use CyteType's hosted models,** also for runs annotated with your own models.
- **Your data still goes through CyteType.** The files are uploaded and the report is hosted as usual (see [Security and privacy](../support/security-and-privacy.md)).

## Configure a model

Pass a list of model configurations: dictionaries in Python, named lists in R.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
import os

adata = annotator.run(
    study_context="...",
    llm_configs=[
        {
            "provider": "openai",
            "name": "<model-name>",
            "apiKey": os.environ["OPENAI_API_KEY"],
        }
    ],
)
```
:::
:::{tab-item} R
:sync: r
```r
pbmc <- CyteTypeR(
  pbmc, prepped_data,
  study_context = "...",
  llm_configs = list(
    list(
      provider = "openai",
      name = "<model-name>",
      apiKey = Sys.getenv("OPENAI_API_KEY")
    )
  ),
  save_query = FALSE  # query.json would otherwise contain the key
)
```

CyteTypeR also accepts a single configuration as a named list, without the outer `list()`.
:::
::::

`OPENAI_API_KEY` stands for a variable you set yourself; CyteType does not read it automatically. Replace `<model-name>` with the model name as your provider spells it.

## Configuration fields

| Field | Required | Meaning |
|---|---|---|
| `provider` | Yes | The provider; see [Providers](#llm-configuration-providers). |
| `name` | Yes | The model name, as your provider spells it. |
| `apiKey` | Yes, except with `bedrock` | Your API key for the provider. |
| `baseUrl` | No | The address of an OpenAI-compatible endpoint. Used by `openai` (and by `vertex`); the other providers connect to their standard endpoint and ignore it. |
| `modelSettings` | No | Settings sent with each request to the model, such as `temperature` or `max_tokens`. |
| `targetAgents` | No | The agents this model serves; see [Several models](#llm-configuration-agents). |
| `awsAccessKeyId`, `awsSecretAccessKey`, `awsDefaultRegion` | With `bedrock`, all three | Your AWS credentials and region; see [AWS Bedrock](#llm-configuration-bedrock). |
| `allowFallback` | No | `True` makes this a fallback model, used only when the other models for the same agent fail. Default `False`. |

The client also accepts `skipValidation`, but the server ignores it and tests every model (see [Validation](#llm-configuration-validation)).

(llm-configuration-providers)=
## Providers

| Provider | `provider` value |
|---|---|
| Anthropic | `anthropic` |
| AWS Bedrock | `bedrock` |
| DeepSeek | `deepseek` |
| Fireworks AI | `fireworks` |
| Google | `google` |
| Groq | `groq` |
| Hugging Face | `huggingface` |
| Mistral | `mistral` |
| Nebius | `nebius` |
| OpenAI, or any OpenAI-compatible endpoint | `openai` |
| OpenRouter | `openrouter` |
| Together AI | `together` |
| Google Vertex AI | `vertex` |
| xAI | `xai` |

For `vertex`, set `baseUrl` to the OpenAI-compatible endpoint of your Vertex AI project; CyteType sends `apiKey` with each request as the `key` parameter.

(llm-configuration-agents)=
## Several models and `targetAgents`

A job is carried out by several agents, each responsible for one part of the work:

| Agent | Role |
|---|---|
| `contextualizer` | Interprets the study context and builds each cluster's context. |
| `annotator` | Proposes and weighs candidate identities, places them in the reference hierarchy and refines labels across clusters. |
| `reviewer` | Reviews each annotation independently: confidence, heterogeneity and the comparison with your own labels. |

`targetAgents` assigns models to agents. It says which agents a model is for; it does not prevent the model from being used elsewhere:

- **A model with `targetAgents`** serves the agents in its list.
- **A model without `targetAgents`** serves every agent.
- **Agents that no model lists** are added to the first model in `llm_configs`. So one model alone covers every agent.
- **Several models for one agent.** CyteType spreads the calls across them and tries another one when a call fails. Models with `allowFallback` set are tried after the others, except models that have just failed several times in a row, which are tried last.

With a single model, the same model annotates and reviews each cluster. To have the review done by a different model, give the `reviewer` its own model.

`targetAgents` also accepts `summarizer`, `clinician` and `chat`, but they have no effect: no step of the current workflow uses the summarizer, and target analysis and Chat always run on CyteType's hosted models.

To dedicate each model to its own agents, give every model a `targetAgents` list:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
llm_configs = [
    {
        "provider": "anthropic",
        "name": "<model-name>",
        "apiKey": os.environ["ANTHROPIC_API_KEY"],
        "targetAgents": ["annotator", "reviewer"],
    },
    {
        "provider": "openai",
        "name": "<model-name>",
        "apiKey": os.environ["OPENAI_API_KEY"],
        "targetAgents": ["contextualizer"],
    },
    {
        "provider": "openrouter",
        "name": "<model-name>",
        "apiKey": os.environ["OPENROUTER_API_KEY"],
        "targetAgents": ["annotator", "reviewer"],
        "allowFallback": True,
    },
]
```
:::
:::{tab-item} R
:sync: r
```r
llm_configs <- list(
  list(
    provider = "anthropic",
    name = "<model-name>",
    apiKey = Sys.getenv("ANTHROPIC_API_KEY"),
    targetAgents = c("annotator", "reviewer")
  ),
  list(
    provider = "openai",
    name = "<model-name>",
    apiKey = Sys.getenv("OPENAI_API_KEY"),
    targetAgents = c("contextualizer")
  ),
  list(
    provider = "openrouter",
    name = "<model-name>",
    apiKey = Sys.getenv("OPENROUTER_API_KEY"),
    targetAgents = c("annotator", "reviewer"),
    allowFallback = TRUE
  )
)
```
:::
::::

Here the first model serves the annotator and the reviewer, with the third as their fallback, and the second model serves the contextualizer.

(llm-configuration-bedrock)=
## AWS Bedrock

Bedrock uses AWS credentials instead of `apiKey`:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
llm_configs = [
    {
        "provider": "bedrock",
        "name": "<bedrock-model-id>",
        "awsAccessKeyId": os.environ["AWS_ACCESS_KEY_ID"],
        "awsSecretAccessKey": os.environ["AWS_SECRET_ACCESS_KEY"],
        "awsDefaultRegion": "<aws-region>",
    }
]
```
:::
:::{tab-item} R
:sync: r
```r
llm_configs <- list(
  list(
    provider = "bedrock",
    name = "<bedrock-model-id>",
    awsAccessKeyId = Sys.getenv("AWS_ACCESS_KEY_ID"),
    awsSecretAccessKey = Sys.getenv("AWS_SECRET_ACCESS_KEY"),
    awsDefaultRegion = "<aws-region>"
  )
)
```
:::
::::

- **`name`** is the model ID as Bedrock lists it.
- **All three AWS fields are required.** If you set one of them, the client asks for the other two.
- **There is no field for a session token,** so temporary AWS credentials cannot be used. Use an access key whose permissions allow calling the model in that region.

(llm-configuration-openai-compatible)=
## OpenAI-compatible endpoints

Many serving tools and providers offer an OpenAI-compatible API. To use one, set `provider` to `openai` and `baseUrl` to the endpoint:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
llm_configs = [
    {
        "provider": "openai",
        "name": "<model-name>",
        "apiKey": os.environ["MY_ENDPOINT_API_KEY"],
        "baseUrl": "https://llm.example.org/v1",
    }
]
```
:::
:::{tab-item} R
:sync: r
```r
llm_configs <- list(
  list(
    provider = "openai",
    name = "<model-name>",
    apiKey = Sys.getenv("MY_ENDPOINT_API_KEY"),
    baseUrl = "https://llm.example.org/v1"
  )
)
```
:::
::::

- **`baseUrl`** is the address of the API, usually ending in `/v1`. With it, CyteType uses the chat completions API. Without it, `openai` connects to OpenAI's own API, so leave it out for OpenAI itself.
- **Reachable from the internet.** CyteType's servers call the endpoint, not your computer, so an address that only works inside your network or on `localhost` does not work.
- **Tool calling.** The model must support tool (function) calling through the endpoint.
- **`apiKey` is required.** If the endpoint does not check keys, pass any non-empty text.

(llm-configuration-ollama)=
## Local models with Ollama

[Ollama](https://ollama.com) serves models on your own machine through an OpenAI-compatible API. Because CyteType's servers call the model, it needs a public address. This tutorial gives a local Ollama server a temporary public address with [ngrok](https://ngrok.com). ngrok is quick to set up and needs no domain name; on its free plan, the public address changes every time you restart ngrok, and the tunnel closes when you close its terminal.

You need:

- a Linux machine with Ollama installed;
- at least one model, downloaded with `ollama pull <model-name>`;
- a model that supports tools. The [Ollama library](https://ollama.com/library) shows which models do;
- enough memory and compute for the model. If the model is too large for your machine, Ollama can answer with `Error 500: Internal Server Error`.

### Step 1: Make Ollama listen on all interfaces

By default, Ollama only accepts requests from the same machine (`127.0.0.1`) on port `11434`. To let a reverse proxy connect to it, configure it to listen on all network interfaces (`0.0.0.0`).

Create an override directory for the Ollama service:

```bash
sudo mkdir -p /etc/systemd/system/ollama.service.d
```

Create a configuration file:

```bash
sudo nano /etc/systemd/system/ollama.service.d/override.conf
```

Add the following content:

```text
[Service]
Environment="OLLAMA_HOST=0.0.0.0"
```

Save the file and exit the editor (in nano, press `Ctrl+O`, then `Enter` to write, then `Ctrl+X` to exit). Reload systemd and restart Ollama:

```bash
sudo systemctl daemon-reload
sudo systemctl restart ollama
```

Check that Ollama listens on all interfaces. Look for an entry with `0.0.0.0:11434` or `*:11434`:

```bash
ss -tlnp | grep 11434
```

The output looks similar to this:

```text
LISTEN 0  4096       0.0.0.0:11434       0.0.0.0:*  users:(("ollama",pid=1234,fd=7))
```

Ollama now accepts connections from other devices on your local network and from the tunnel.

:::{note}
On WSL without systemd, start Ollama in a separate terminal instead:

```bash
OLLAMA_HOST=0.0.0.0 ollama serve
```

Ollama keeps a model in memory for 5 minutes after a response. To keep it loaded for the whole job, start it with `OLLAMA_HOST=0.0.0.0 OLLAMA_KEEP_ALIVE=-1 ollama serve`.
:::

### Step 2: Install ngrok

On Debian or Ubuntu:

```bash
curl -sSL https://ngrok-agent.s3.amazonaws.com/ngrok.asc \
  | sudo tee /etc/apt/trusted.gpg.d/ngrok.asc >/dev/null \
  && echo "deb https://ngrok-agent.s3.amazonaws.com buster main" \
  | sudo tee /etc/apt/sources.list.d/ngrok.list \
  && sudo apt update \
  && sudo apt install ngrok
```

For other systems, see the [ngrok download page](https://ngrok.com/downloads/).

### Step 3: Sign up and add your authtoken

Create a free account on the [ngrok dashboard](https://dashboard.ngrok.com/), copy your authtoken from there, and add it to your local configuration (once):

```bash
ngrok config add-authtoken <token>
```

### Step 4: Start the tunnel

Forward a public address to Ollama's port:

```bash
ngrok http 11434
```

ngrok shows your public address on the `Forwarding` line:

```text
Session Status                online
Account                       Your Name (Plan: Free)
Version                       3.x.x
Region                        United States (us)
Forwarding                    https://random-string-of-letters.ngrok-free.app -> http://localhost:11434
```

Keep this terminal open; closing it stops the tunnel.

### Step 5: Test the endpoint

The OpenAI-compatible API is under `/v1`. Replace the address with yours and `<model-name>` with your model:

```bash
curl https://your-ngrok-url.ngrok-free.app/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "<model-name>",
    "messages": [
      {
        "role": "user",
        "content": "Why is the sky blue?"
      }
    ]
  }'
```

### Step 6: Use the model in CyteType

Use the `openai` provider with your public address as `baseUrl`. Ollama does not check API keys, so any text works as `apiKey`:

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(
    study_context="...",
    llm_configs=[
        {
            "provider": "openai",
            "name": "<model-name>",
            "apiKey": "ollama",
            "baseUrl": "https://your-ngrok-url.ngrok-free.app/v1",
        }
    ],
)
```
:::
:::{tab-item} R
:sync: r
```r
pbmc <- CyteTypeR(
  pbmc, prepped_data,
  study_context = "...",
  llm_configs = list(
    list(
      provider = "openai",
      name = "<model-name>",
      apiKey = "ollama",
      baseUrl = "https://your-ngrok-url.ngrok-free.app/v1"
    )
  )
)
```
:::
::::

Before you submit:

- **Keep Ollama and the tunnel running until the job completes.** CyteType tests the model when you submit and calls it throughout the job.
- **Keep `n_parallel_clusters` low,** at the default of 2 or at 1, because a single local server answers few requests at a time.

:::{warning}
Anyone who knows the ngrok address can send requests to your model while the tunnel is open. Stop ngrok when the job is done.
:::

(llm-configuration-validation)=
## Validation before the job starts

Your configuration is checked twice before the job starts:

1. **By the client, before anything is uploaded.** It checks the fields: a known provider, a model name, and either `apiKey` or all three AWS fields, and valid agent names in `targetAgents`. Python raises a `pydantic.ValidationError`; CyteTypeR stops with a message such as `provider must be one of: ...`.
2. **By the server, at submission.** It sends each model a short test that needs a tool call, with up to three attempts. If any model fails, the whole request is rejected: Python raises `LLMValidationError`, and CyteTypeR stops with `Job submission failed: [LLM_VALIDATION_FAILED] ...`. The message names the models that failed.

A rejected request creates no job and does not count toward your limits. The files had already been uploaded, though, and are uploaded again when you resubmit. The usual causes are a misspelled model name, a wrong or expired key, an endpoint that cannot be reached from the internet, and a model that does not support tool calling.

(llm-configuration-keys)=
## Keep keys out of files and code

- **Read keys from environment variables or a secret store** rather than writing them into scripts or notebooks, as in the examples on this page.
- **Python** never writes `llm_configs` or their keys to `query.json`, and does not store them in the AnnData.
- **CyteTypeR** writes `llm_configs`, including their keys, to `query.json` when `save_query = TRUE`, the default. Pass `save_query = FALSE` whenever you use your own models, and delete any `query.json` that contains keys.
- **Your keys are sent to CyteType** with each job, so that its servers can call your provider. Consider a separate key for CyteType that you can revoke on its own.

(llm-configuration-parallel)=
## Parallel requests

`n_parallel_clusters` (from 1 to 50, default 2) sets how many clusters are annotated at the same time when you use your own models. Higher values finish sooner but send more requests at once, which can run into your provider's rate limits. Without your own models, the server uses its own setting and ignores this value.

::::{tab-set}
:::{tab-item} Python
:sync: python
```python
adata = annotator.run(study_context="...", llm_configs=llm_configs, n_parallel_clusters=8)
```
:::
:::{tab-item} R
:sync: r
```r
pbmc <- CyteTypeR(
  pbmc, prepped_data,
  study_context = "...",
  llm_configs = llm_configs,
  n_parallel_clusters = 8L,
  save_query = FALSE
)
```
:::
::::

## See also

- [LLM configuration reference](../reference/api/llm-config.html): every field of a configuration, as the Python client validates it.
- [Usage limits and quotas](../platform/limits.md): the limits that apply to every job, including jobs with your own models.
- [Security and privacy](../support/security-and-privacy.md): what passes through CyteType when you use your own provider.
